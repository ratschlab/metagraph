#!/usr/bin/env python3
"""Local library runner of area G (DESIGN-v8 §2, §7 item 2; FORMAT-v8 §1-§3): the traversal library at a given
commit, run over the run's graphlet documents with the manifest's op plans, one isolated process per way, writing
FORMAT §3's per-call records as a runner writes its executions.

    run_graphlet_local.py --run-dir RUN --manifest M --library-root DIR [--ways W1,W2,W4] [--scenarios ID,...]
                          [--limits L.json] [--out RUN] [--repeat 0] [--mgt-template T]
                          [--ways drift --drift-library-root DIR2]
    run_graphlet_local.py --doc [ID=]PATH [--doc ...] --plan PLAN.json --library-root DIR --out DIR [--ways W4]
    run_graphlet_local.py --export-library COMMIT --to DIR [--repo REPO]
    run_graphlet_local.py --drift-report --out RUN

Ways and their documents (DESIGN-v8 §2):
  W1 -> arms/dl-lib     the direct-local answer arms/direct-local/raw/<id>.direct-local.<repeat>.json.gz (answer.body:
                        the /traverse response), turned into the service's document form by the library's own
                        parser.standalone_text(results[0].graphlet, results[0], response) in the way's process (byte
                        for byte dump(from_response(...), envelope=True), parser.py:1343-1366)
  W2 -> arms/dr-lib     the same from arms/direct-remote
  W4 -> arms/mcp-lib    the MCP's exported file arms/mcp/docs/<id>.mgt.gz, as it is (FORMAT v2 C1: one MCP execution
                        per area G scenario, no repeat; --mgt-template overrides it)
  drift -> arms/drift-lib  W4's documents and ops under --drift-library-root (the branch head); report-only, listed by
                        --drift-report against W4's records
Offline (--doc/--plan: the tests, the rehearsal): a .json/.json.gz document (a /traverse response or a raw record)
enters as W1's does, any other (.mgt, .mgt.gz, a bare server body too) as W4's; every way named runs every document.

Every document enters its way the same way: written to the way's spool and loaded by the library's graphlet_load (the
load path of a standalone file, the one the service names for its links). The way's GraphletTools is configured from
--limits (default max_bytes 2048, sequence_max_bytes 16384, local_limits null, max_graphlet_mb 1e6: DESIGN §3), with
a fixed cursor secret, and the limits are recorded with every call.

Process isolation: one child process per way, `python -s -B [-P] run_graphlet_local.py --child`, its environment
scrubbed to PATH/HOME/TMPDIR/LANG/LC_* plus the way's PYTHONHASHSEED (distinct per way, recorded), its sys.path the
library root first (the script's directory and the working directory removed); the child imports nothing of the
suite. Not `python3 -I`: -I implies -E, which ignores PYTHONHASHSEED (measured on CPython 3.15), and the per-way hash
seed is what makes a set or dict order leaking into an answer show up as a W1/W2/W4 difference (DESIGN §2, D-R2).
The child checks that the library it imported lies under the root, and the parent that the seed took effect.

Supervision (VERIFY-L VL3): every call waits for its reply with select() up to the op's `timeout_s` (else
--op-timeout, default 120 s). A call that does not answer in time is recorded `timeout` {timeout_s, page}, a process
that ends during a call `child_died` {exit, page}; either way the process is killed and restarted, the documents
loaded so far reloaded into it, and the next op runs (handles an op made are gone with the old process). An answer
that is no JSON is that call's `exception` {type}, one holding a NaN or an infinity `non_finite` {pointer, page}
(VL9: the records are strict JSON, as the service refuses such numbers with number_not_finite). The child reads its
requests from a duplicate of fd 0 and points fd 0 at /dev/null, as it points fd 1 at its log. A way whose process does not start, or restarts more than 20
times, is reported and the next way runs; the invocation then exits 1.

Op execution (FORMAT §1-§2): per way, scenarios in manifest order, each scenario's ops in order; an op runs on a way
when its `ways` (else the scenario's graphlet.ways.default, else every way) names the way (drift runs W4's ops).
Arguments are resolved by graphlet_select.resolve from this way's earlier assembled answers; {"graphlet": "<id>"}
anywhere in them is that scenario's graphlet, loaded in the same process (compare's `b`, export's `other`); the
graphlet's own handle is filled in as `handle` (compare: `a`) when the plan names none. Paging and assembly are
graphlet_paging's (FORMAT v2 A2, shared with async's runner): "none" one call, "all" every page (next_cursor;
graphlet_sequence: next_from as start), {"pages": N} at most N, the pages assembled into one answer.

Handles in records: a way's random `g_` handles are written as symbols, relative to the op's scenario: "<handle>" the
scenario's own graphlet (FORMAT's args_resolved form), "<handle:S>" scenario S's, "<handle:OP>" the one op OP of this
scenario made (graphlet_subtrie, graphlet_load), "<handle:S/OP>" one another scenario's op made, "<handle:?>" any
other; cursors in args_resolved as "<cursor>"; the export and spool directories as "<export_dir>", "<spool_dir>".
The raw pages keep everything as the library answered it.

Files of a way's directory: graphlet_ops.jsonl (one record per op, FORMAT §3, plus `paging_complete` and, for
graphlet_export / graphlet_save, `file`), graphlet_docs.jsonl (one line per document: its source, digests, load
receipt or error), raw/<scenario>.<op>.<page>.json.gz (page 0 first) and raw/<scenario>.<op>.file<.ext>.gz (an
export's text, handles as symbols), docs/<scenario>.mgt.gz (the document the way processed), run.json, child.log,
work/ (the spool and the export directory). An existing graphlet_ops.jsonl is never appended to: use another --out.
"""

import argparse
import gzip
import hashlib
import json
import os
import platform
import re
import select
import subprocess
import sys
import time
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
CHILD = "--child"

WAYS = ("W1", "W2", "W4", "drift")
WAY_DIRS = {"W1": "dl-lib", "W2": "dr-lib", "W4": "mcp-lib", "drift": "drift-lib"}
WAY_SOURCE = {"W1": "direct-local", "W2": "direct-remote", "W4": "mcp", "drift": "mcp"}
# the plan's way an op must name to run on this way (the drift lane re-runs W4's ops at another commit)
PLAN_WAY = {"W1": "W1", "W2": "W2", "W4": "W4", "drift": "W4"}
DEFAULT_HASHSEEDS = {"W1": 1001, "W2": 2002, "W4": 4004, "drift": 5005}
# FORMAT v2 C1: the MCP's exported document of a scenario, which W4 and drift read; each way writes the document it
# processed to arms/<way dir>/docs/<scenario>.mgt.gz
DEFAULT_MGT_TEMPLATE = "arms/mcp/docs/{id}.mgt.gz"
DEFAULT_LIMITS = {"max_bytes": 2048, "sequence_max_bytes": 16384, "local_limits": None, "max_graphlet_mb": 1e6}
# fixed, so that a cursor (raw pages only) is the same string in every process of one commit
CURSOR_SECRET = "metagraph-suite-v8-local-cursors"
ENV_KEEP = ("PATH", "HOME", "TMPDIR", "LANG", "LC_ALL", "LC_CTYPE")
# supervision of a way's process (VERIFY-L VL3): a call that answers within no time limit is a `timeout`, a process
# that ends in a call `child_died`; either restarts the process (its documents reloaded) and the run goes on. An op's
# own limit is its manifest field `timeout_s` (the rehearsal's per-op medians, DESIGN-v8 §4), else --op-timeout; the
# heavy lane's 120 s is the default (DESIGN-v8 §3: every op <= 1/3 of its lane)
DEFAULT_OP_TIMEOUT_S = 120.0
STARTUP_TIMEOUT_S = 120.0
MAX_RESTARTS = 20

# the argument naming the graphlet a tool reads; graphlet_list, graphlet_load and the backend tools read none
HANDLE_ARG = dict([(t, "handle") for t in (
    "graphlet_summary", "graphlet_walks", "graphlet_walk", "graphlet_support", "graphlet_labels", "graphlet_splits",
    "graphlet_claims", "graphlet_sequence", "graphlet_export", "graphlet_subtrie", "graphlet_save", "graphlet_free",
    "traverse_continue")] + [("graphlet_compare", "a")])
# the tools whose answer's `handle` is a new one (subtrie :2390, load :2470; fetch and continue need a backend)
MINTING = ("graphlet_subtrie", "graphlet_load", "traverse_fetch", "traverse_continue")
FILE_TOOLS = ("graphlet_export", "graphlet_save")
HANDLE_RE = re.compile(r"g_[0-9a-f]{12}")
SYMBOL_RE = re.compile(r"<handle(?::([^<>]+))?>\Z")


# ================================================================================================ the way's process
def child_main(argv):
    """The process of one way: reads one JSON request per line on stdin, answers one JSON line on stdout. Requests:
    init {limits, spool_dir, export_dir, secret} -> the interpreter, hash seed and library it runs; convert {response,
    index} -> {text} (parser.standalone_text) or {no_graphlet}; call {tool, args} -> {answer, call_ms} or
    {unknown_tool}; quit. An exception (a library bug) is answered {ok: false, type, error, traceback}, never fatal."""
    root = os.path.realpath(argv[argv.index("--library-root") + 1])
    here = os.path.realpath(HERE)
    sys.path[:] = [p for p in sys.path if p and os.path.realpath(p) != here]
    sys.path.insert(0, root)
    channel = os.fdopen(os.dup(1), "wb")
    requests = os.fdopen(os.dup(0), "rb")
    os.dup2(2, 1)          # whatever the library prints goes to the log, never into the protocol
    sys.stdout = sys.stderr
    null = os.open(os.devnull, os.O_RDONLY)
    os.dup2(null, 0)       # and whatever it reads from stdin is nothing, never a request
    os.close(null)
    sys.stdin = open(os.devnull, "r")
    state = {"root": root}
    for line in requests:
        req = json.loads(line.decode("utf-8"))
        try:   # the answer serialised inside: one that is no JSON is that call's exception, not the way's end
            result = _child_do(state, req)
            where = _non_finite(result.get("answer") if isinstance(result, dict) else None)
            if where is not None:   # NaN / Infinity: no JSON a strict reader takes (VERIFY-L VL9)
                data = json.dumps({"ok": False, "non_finite": where})
            else:
                data = json.dumps({"ok": True, "result": result}, ensure_ascii=False, allow_nan=False)
        except Exception as e:  # a library exception is the op's error (code exception), not the way's end
            data = json.dumps({"ok": False, "type": type(e).__name__, "error": str(e),
                               "traceback": traceback.format_exc(limit=12)}, ensure_ascii=False)
        channel.write((data + "\n").encode("utf-8"))
        channel.flush()
        if req.get("cmd") == "quit":
            break
    return 0


def _non_finite(obj, path=""):
    """The RFC 6901 pointer of the first NaN or infinite number in |obj|, or None."""
    if isinstance(obj, float) and (obj != obj or obj in (float("inf"), float("-inf"))):
        return path
    items = obj.items() if isinstance(obj, dict) else enumerate(obj) if isinstance(obj, list) else ()
    for k, v in items:
        got = _non_finite(v, "%s/%s" % (path, str(k).replace("~", "~0").replace("/", "~1")))
        if got is not None:
            return got
    return None


def _tool_limits(spec):
    """The limits file's local_limits -> ToolLimits (None: the tools run unbudgeted; {} or true: the defaults; else
    {view|heavy|parse|*_ceiling: {work_units, memory_mb, deadline_s}, per_tool: {tool: {...}}})."""
    if spec is None or spec is False:
        return None
    from metagraph.traverse.budget import LocalLimits
    from metagraph.traverse.mcp_tools import ToolLimits
    kw = {}
    for key, value in ({} if spec is True else spec).items():
        if key == "per_tool":
            kw[key] = dict((t, LocalLimits(**v)) for t, v in value.items())
        else:
            kw[key] = LocalLimits(**value)
    return ToolLimits(**kw)


def _child_do(state, req):
    cmd = req.get("cmd")
    if cmd == "init":
        import metagraph.traverse as T
        from metagraph.traverse import budget, parser
        from metagraph.traverse.mcp_tools import GraphletTools, tool_names
        lib = os.path.realpath(os.path.dirname(T.__file__))
        if not lib.startswith(state["root"] + os.sep):
            raise RuntimeError("metagraph.traverse was imported from %s, not from the library root %s"
                               % (lib, state["root"]))
        limits = req["limits"]
        kw = {"secret": req["secret"].encode("utf-8"), "export_dir": req["export_dir"],
              "max_bytes": limits["max_bytes"], "sequence_max_bytes": limits["sequence_max_bytes"]}
        local = _tool_limits(limits.get("local_limits"))
        if local is not None:
            kw["local_limits"] = local
        state["tools"] = GraphletTools(T.GraphletStore(req["spool_dir"]), **kw)
        state["names"] = set(tool_names())
        state["parser"] = parser
        return {"python": platform.python_version(), "implementation": platform.python_implementation(),
                "executable": sys.executable, "hashseed": os.environ.get("PYTHONHASHSEED"),
                "hash_randomization": sys.flags.hash_randomization, "hash_probe": hash("metagraph suite v8"),
                "work_model": getattr(budget, "WORK_MODEL", None), "library_file": lib,
                "tools": sorted(state["names"])}
    if cmd == "convert":
        resp, i = req["response"], req.get("index", 0)
        results = resp.get("results") if isinstance(resp, dict) else None
        if not isinstance(results, list) or i >= len(results) or not isinstance(results[i], dict):
            return {"refused": "no_results", "text": "the response has no results[%d]" % i}
        r = results[i]
        if not isinstance(r.get("graphlet"), str):
            if "error" in r:   # a failed seed (SPEC-core §7.1's failed shape {seed, error}): the one no_graphlet
                return {"no_graphlet": "results[%d] failed: %s" % (i, r["error"])}
            return {"refused": "no_graphlet_no_error",
                    "text": "results[%d] carries neither a graphlet nor an error" % i}
        return {"text": state["parser"].standalone_text(r["graphlet"], r, resp)}
    if cmd == "call":
        if req["tool"] not in state["names"]:
            return {"unknown_tool": req["tool"]}
        fn = getattr(state["tools"], req["tool"])
        t0 = time.perf_counter()
        answer = fn(**req["args"])
        return {"answer": answer, "call_ms": round((time.perf_counter() - t0) * 1000.0, 3)}
    if cmd == "quit":
        return None
    raise ValueError("unknown command %r" % (cmd,))


if sys.argv[1:2] != [CHILD]:   # the runner (and the tests): the suite's own modules; the child imports none of them
    sys.path.insert(0, HERE)
    import graphlet_doc  # noqa: E402
    import graphlet_paging as gp  # noqa: E402
    import graphlet_select  # noqa: E402
    import suite_common as sc  # noqa: E402


# ================================================================================================ helpers
def safe(name):
    """A scenario or op id as a file name part (suite_common.raw_name's rule)."""
    return re.sub(r"[^A-Za-z0-9._-]", "_", str(name))


def iso(ts):
    return sc.iso_utc(ts)


def git(repo, *args):
    return subprocess.check_output(("git", "-C", repo) + args, stderr=subprocess.DEVNULL).decode("utf-8").strip()


def library_path(root):
    """The sys.path entry of a library root: DIR/metagraph/api/python (a checkout or an --export-library copy), or
    DIR itself when it holds metagraph/traverse."""
    root = os.path.abspath(root)
    for cand in (os.path.join(root, "metagraph", "api", "python"), root):
        if os.path.isfile(os.path.join(cand, "metagraph", "traverse", "mcp_tools.py")):
            return cand
    raise SystemExit("run_graphlet_local: %s holds no metagraph/api/python/metagraph/traverse" % root)


def library_commit(root, explicit=None):
    """The commit of a library root: --library-commit, else the LIBRARY_COMMIT file --export-library writes, else the
    checkout's HEAD (with "+dirty" when metagraph/api/python differs from it), else None."""
    if explicit:
        return explicit
    path = library_path(root)
    for d in (os.path.abspath(root), os.path.normpath(os.path.join(path, "..", "..", ".."))):
        marker = os.path.join(d, "LIBRARY_COMMIT")
        if os.path.isfile(marker):
            with open(marker, "r", encoding="utf-8") as fh:
                return fh.read().strip()
    try:
        top = git(path, "rev-parse", "--show-toplevel")
        head = git(top, "rev-parse", "HEAD")
        dirty = git(top, "status", "--porcelain", "--", os.path.relpath(path, top))
    except (OSError, subprocess.CalledProcessError):
        return None
    return head + ("+dirty" if dirty else "")


def python_version(python):
    out = subprocess.check_output([python, "-c", "import sys; print('%d %d' % sys.version_info[:2])"],
                                  env=dict((k, os.environ[k]) for k in ENV_KEEP if k in os.environ))
    return tuple(int(x) for x in out.split())


def read_limits(path):
    """DEFAULT_LIMITS updated by the limits file (an unknown key is refused)."""
    limits = dict(DEFAULT_LIMITS)
    if path:
        given = sc.read_json(path)
        if not isinstance(given, dict) or set(given) - set(DEFAULT_LIMITS):
            raise SystemExit("run_graphlet_local: %s is not a limits object (keys %s)"
                             % (path, ", ".join(sorted(DEFAULT_LIMITS))))
        limits.update(given)
    return limits


def read_text(path):
    """A document's text: gunzipped when the file is gzip, UTF-8 strictly."""
    with open(path, "rb") as fh:
        data = fh.read()
    if data[:2] == b"\x1f\x8b":
        data = gzip.decompress(data)
    return data.decode("utf-8")


# ============================================================================================ the child, parent side
class ChildTimeout(Exception):
    """The way's process gave no answer within the call's time limit."""


class ChildDied(Exception):
    """The way's process ended during a call."""

    def __init__(self, exit_code):
        super().__init__("the way's process ended (exit %s)" % exit_code)
        self.exit = exit_code


class Supervised(Exception):
    """A call the supervision ended (timeout, child_died): the process was restarted, the call has no answer."""

    def __init__(self, code, fields):
        super().__init__(code)
        self.code = code
        self.fields = fields


class WayFailed(Exception):
    """A way that cannot go on (its process does not start, or restarts too often): later ways still run."""


class Child(object):
    """One way's process (child_main), driven over its stdin/stdout; a reply is awaited with select() up to the
    call's time limit (ChildTimeout), a closed stdout is ChildDied."""

    def __init__(self, python, root, hashseed, cwd, log_path):
        flags = ["-s", "-B"] + (["-P"] if python_version(python) >= (3, 11) else [])
        self.cmd = [python] + flags + [os.path.abspath(__file__), CHILD, "--library-root", root]
        self.env = dict((k, os.environ[k]) for k in ENV_KEEP if k in os.environ)
        self.env["PYTHONHASHSEED"] = str(hashseed)
        self.log = open(log_path, "ab")
        self.proc = subprocess.Popen(self.cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=self.log,
                                     cwd=cwd, env=self.env, bufsize=0)
        self.buf = b""

    def ask(self, req, timeout=None):
        try:
            self.proc.stdin.write((json.dumps(req, ensure_ascii=False) + "\n").encode("utf-8"))
            self.proc.stdin.flush()
        except (BrokenPipeError, OSError):
            raise ChildDied(self.proc.wait()) from None
        fd = self.proc.stdout.fileno()
        deadline = None if timeout is None else time.monotonic() + timeout
        while b"\n" not in self.buf:
            wait = None if deadline is None else deadline - time.monotonic()
            if wait is not None and wait <= 0:
                raise ChildTimeout()
            ready, _, _ = select.select([fd], [], [], wait)
            if not ready:
                raise ChildTimeout()
            chunk = os.read(fd, 1 << 20)
            if not chunk:
                raise ChildDied(self.proc.wait())
            self.buf += chunk
        line, _, self.buf = self.buf.partition(b"\n")
        return json.loads(line.decode("utf-8"))

    def close(self, kill=False):
        if self.proc.poll() is None and not kill:
            try:
                self.ask({"cmd": "quit"}, 10.0)
            except (ChildTimeout, ChildDied, OSError, ValueError):
                pass
            self.proc.stdin.close()
            try:
                self.proc.wait(timeout=30)
            except subprocess.TimeoutExpired:
                kill = True
        if kill and self.proc.poll() is None:
            self.proc.kill()
            self.proc.wait()
        for pipe in (self.proc.stdin, self.proc.stdout):
            if not pipe.closed:
                pipe.close()
        self.log.close()


# ================================================================================================ one way
class Way(object):
    def __init__(self, runner, name):
        self.r = runner
        self.name = name
        self.dir = os.path.join(runner.out, "arms", WAY_DIRS[name])
        self.ops_path = os.path.join(self.dir, "graphlet_ops.jsonl")
        if os.path.exists(self.ops_path):
            raise SystemExit("run_graphlet_local: %s exists; records are never appended to: use another --out"
                             % self.ops_path)
        for sub in ("raw", "docs", os.path.join("work", "spool", "inputs"), os.path.join("work", "exports")):
            os.makedirs(os.path.join(self.dir, sub), exist_ok=True)
        self.root = runner.roots[name]
        self.hashseed = runner.hashseeds[name]
        self.spool = os.path.realpath(os.path.join(self.dir, "work", "spool"))
        self.exports = os.path.realpath(os.path.join(self.dir, "work", "exports"))
        self.child = None
        self.restarts = []   # [{reason, utc}] of every restart of the way's process (supervision)
        self.start_child()
        self.library = {"commit": runner.commits[name], "python": self.info["python"], "hashseed": self.hashseed,
                        "work_model": self.info["work_model"], "limits": runner.limits}
        self.handles = {}    # scenario -> its graphlet's handle in this process
        self.minted = {}     # (scenario, op) -> the handle that op made
        self.owner = {}      # handle -> (scenario, op or None)
        self.docs = {}       # scenario -> its graphlet_docs line
        self.counts = {"ops": 0, "answered": 0, "error": 0, "selector_error": 0, "calls": 0}
        self.started = time.time()
        self.ops_file = open(self.ops_path, "w", encoding="utf-8")
        self.docs_file = open(os.path.join(self.dir, "graphlet_docs.jsonl"), "w", encoding="utf-8")
        self.write_run("running")

    # -- the way's process ------------------------------------------------------------------------------------
    def start_child(self):
        """Start (or restart) the way's process and initialise the library in it; WayFailed when it does not start."""
        self.child = Child(self.r.python, library_path(self.root), self.hashseed, os.path.join(self.dir, "work"),
                           os.path.join(self.dir, "child.log"))
        try:
            msg = self.child.ask({"cmd": "init", "limits": self.r.limits, "spool_dir": self.spool,
                                  "export_dir": self.exports, "secret": CURSOR_SECRET}, STARTUP_TIMEOUT_S)
        except (ChildTimeout, ChildDied) as e:
            self.child.close(kill=True)
            raise WayFailed("[%s] the way's process does not start: %s" % (self.name, type(e).__name__)) from None
        if not msg["ok"]:
            self.child.close()
            raise WayFailed("[%s] the library at %s does not start: %s: %s"
                            % (self.name, self.root, msg["type"], msg["error"]))
        info = msg["result"]
        if info["hashseed"] != str(self.hashseed) or info["hash_randomization"] != 1:
            self.child.close()
            raise WayFailed("[%s] the child runs with hash seed %r (randomization %r), not %d"
                            % (self.name, info["hashseed"], info["hash_randomization"], self.hashseed))
        self.info = info

    def restart(self, reason):
        """After a timeout or a death: a new process, every document loaded before reloaded into it (the handles a
        scenario's ops made are gone: later ops naming them get the library's answer to an unknown handle)."""
        self.restarts.append({"reason": reason, "utc": iso(time.time())})
        self.child.close(kill=True)
        if len(self.restarts) > MAX_RESTARTS:
            raise WayFailed("[%s] the way's process was restarted %d times" % (self.name, len(self.restarts) - 1))
        self.start_child()
        self.minted.clear()
        for sid in list(self.handles):
            path = os.path.join(self.spool, "inputs", safe(sid) + ".mgt")
            try:
                msg = self.child.ask({"cmd": "call", "tool": "graphlet_load", "args": {"path": path}},
                                     self.r.op_timeout)
            except (ChildTimeout, ChildDied) as e:
                raise WayFailed("[%s] reloading %s after a restart: %s" % (self.name, sid, type(e).__name__)) from None
            answer = msg["result"].get("answer") if msg["ok"] else None
            if not isinstance(answer, dict) or not isinstance(answer.get("handle"), str):
                raise WayFailed("[%s] reloading %s after a restart failed" % (self.name, sid))
            self.handles[sid] = answer["handle"]
            self.owner[answer["handle"]] = (sid, None)
        self.write_run("running")

    def ask(self, req, timeout):
        """The child's reply to |req| within |timeout| seconds; a timeout or a death restarts the process and raises
        Supervised(code, fields) for the call in flight (its op is recorded with that error, the run goes on)."""
        try:
            return self.child.ask(req, timeout)
        except ChildTimeout:
            code, fields = "timeout", {"timeout_s": timeout}
        except ChildDied as e:
            code, fields = "child_died", {"exit": e.exit}
        self.restart(code)
        raise Supervised(code, fields)

    # -- records ----------------------------------------------------------------------------------------------
    @staticmethod
    def write_line(fh, obj):
        fh.write(json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False) + "\n")
        fh.flush()

    def write_run(self, state):
        sc.write_json_atomic(os.path.join(self.dir, "run.json"), {
            "run_id": self.r.run_id, "way": self.name, "way_dir": WAY_DIRS[self.name], "state": state,
            "source": self.r.source_note(self.name), "library": dict(self.library, root=os.path.abspath(self.root),
                                                                      file=self.info["library_file"]),
            "interpreter": dict((k, self.info[k]) for k in ("python", "implementation", "executable", "hashseed",
                                                            "hash_randomization", "hash_probe")),
            "child": {"argv": self.child.cmd, "env_keys": sorted(self.child.env), "restarts": self.restarts},
            "limits": self.r.limits, "counts": self.counts, "started_utc": iso(self.started),
            "updated_utc": iso(time.time()), "runner": os.path.basename(__file__)})

    def close(self, state="done"):
        self.child.close(kill=state != "done")
        for fh in (self.ops_file, self.docs_file):
            fh.flush()
            os.fsync(fh.fileno())
            fh.close()
        self.write_run(state)

    # -- handles ----------------------------------------------------------------------------------------------
    def symbol(self, real, scenario):
        own = self.owner.get(real)
        if own is None:
            return "<handle:?>"
        sid, op = own
        if op is None:
            return "<handle>" if sid == scenario else "<handle:%s>" % sid
        return "<handle:%s>" % op if sid == scenario else "<handle:%s/%s>" % (sid, op)

    def symbolise(self, obj, scenario):
        """|obj| with every handle written as its symbol and the way's directories as <export_dir>, <spool_dir>."""
        if isinstance(obj, dict):
            return dict((k, self.symbolise(v, scenario)) for k, v in obj.items())
        if isinstance(obj, list):
            return [self.symbolise(v, scenario) for v in obj]
        if isinstance(obj, str):
            return HANDLE_RE.sub(lambda m: self.symbol(m.group(0), scenario), self.unplace(obj))
        return obj

    def unplace(self, text):
        return text.replace(self.exports + os.sep, "<export_dir>/").replace(self.spool + os.sep, "<spool_dir>/")

    def real(self, value, scenario):
        """A symbol (exactly) -> this process's handle; any other value as it is (the library judges it)."""
        m = SYMBOL_RE.match(value) if isinstance(value, str) else None
        if m is None:
            return value
        name = m.group(1)
        if name is None:
            return self.handles.get(scenario, value)
        if "/" in name:
            return self.minted.get(tuple(name.split("/", 1)), value)
        if (scenario, name) in self.minted:
            return self.minted[(scenario, name)]
        if self.r.has_document(name):
            self.document(name)
            return self.handles.get(name, value)
        return value

    def unsymbolise(self, obj, scenario):
        if isinstance(obj, dict):
            return dict((k, self.unsymbolise(v, scenario)) for k, v in obj.items())
        if isinstance(obj, list):
            return [self.unsymbolise(v, scenario) for v in obj]
        return self.real(obj, scenario)

    def graphlets(self, obj, scenario, where="args"):
        """{"graphlet": "<scenario id>"} at any depth -> that scenario's handle symbol (its document loaded here); a
        scenario without a document, or whose graphlet did not load in this way, is a SelectorError naming |where|."""
        if isinstance(obj, dict) and set(obj) == {"graphlet"}:
            sid = obj["graphlet"]
            if not isinstance(sid, str) or not self.r.has_document(sid):
                raise graphlet_select.SelectorError("%s: {\"graphlet\": %r} names no scenario with a document"
                                                    % (where, sid))
            doc = self.document(sid)
            if sid not in self.handles:
                raise graphlet_select.SelectorError("%s: scenario %s's graphlet is not loaded in this way (%s)"
                                                    % (where, sid, (doc.get("error") or {}).get("code")))
            return self.symbol(self.handles[sid], scenario)
        if isinstance(obj, dict):
            return dict((k, self.graphlets(v, scenario, "%s.%s" % (where, k))) for k, v in obj.items())
        if isinstance(obj, list):
            return [self.graphlets(v, scenario, "%s[%d]" % (where, i)) for i, v in enumerate(obj)]
        return obj

    # -- documents --------------------------------------------------------------------------------------------
    def document(self, sid):
        """Scenario |sid|'s document in this way, loaded once: converted (a response) or read (an MGT file), its
        digests, a copy under docs/, loaded by graphlet_load. -> its graphlet_docs line."""
        if sid in self.docs:
            return self.docs[sid]
        src = self.r.source(self.name, sid)
        line = {"run_id": self.r.run_id, "scenario": sid, "way": self.name, "source": src["rel"], "doc": None,
                "doc_error": None, "load": None, "error": None, "library": self.library}
        self.docs[sid] = line
        try:
            text = self.document_text(src, line)
            if text is not None:
                self.load(sid, text, line)
        finally:
            self.write_line(self.docs_file, line)
        return line

    def document_text(self, src, line):
        if not src["path"] or not os.path.isfile(src["path"]):
            line["error"] = {"code": "no_document", "fields": {"source": src["rel"]}}
            return None
        try:
            if src["kind"] == "mgt":
                return read_text(src["path"])
            obj = json.loads(read_text(src["path"]))
        except (OSError, ValueError) as e:
            line["error"] = {"code": "document_refused", "fields": {"stage": "read", "type": type(e).__name__}}
            line["error_text"] = str(e)
            return None
        resp = obj
        if isinstance(obj, dict) and isinstance(obj.get("answer"), dict):
            # a direct arm's raw record (suite_common.normalise_direct / transport_failure): only a 2xx answer without
            # an error is a response to convert (VERIFY-L VL2); a transport failure (no HTTP status) left no document,
            # an HTTP error or an error body refused it -- neither is a failed seed (no_graphlet)
            answer = obj["answer"]
            http, err = answer.get("http"), answer.get("error")
            code = err.get("code") if isinstance(err, dict) else None
            if not http:
                line["error"] = {"code": "no_document", "fields": {"source": src["rel"], "stage": "fetch",
                                                                   "code": code}}
                line["error_text"] = err.get("text") if isinstance(err, dict) else None
                return None
            if not isinstance(http, int) or not 200 <= http < 300 or err is not None:
                line["error"] = {"code": "document_refused", "fields": {"stage": "fetch", "http": http, "code": code}}
                line["error_text"] = err.get("text") if isinstance(err, dict) else None
                return None
            resp = answer.get("body")
        try:
            msg = self.ask({"cmd": "convert", "response": resp, "index": 0}, self.r.op_timeout)
        except Supervised as e:
            line["error"] = {"code": e.code, "fields": dict(e.fields, stage="convert")}
            return None
        if not msg["ok"]:
            line["error"] = {"code": "document_refused", "fields": {"stage": "convert", "type": msg["type"]}}
            line["error_text"] = msg["error"]
            return None
        if "refused" in msg["result"]:
            line["error"] = {"code": "document_refused", "fields": {"stage": "convert",
                                                                    "reason": msg["result"]["refused"]}}
            line["error_text"] = msg["result"]["text"]
            return None
        if "no_graphlet" in msg["result"]:
            line["error"] = {"code": "no_graphlet", "fields": {}}
            line["error_text"] = msg["result"]["no_graphlet"]
            return None
        return msg["result"]["text"]

    def load(self, sid, text, line):
        data = text.encode("utf-8")
        with gzip.open(os.path.join(self.dir, "docs", safe(sid) + ".mgt.gz"), "wb") as fh:
            fh.write(data)
        try:
            line["doc"] = graphlet_doc.digests(text)
        except ValueError as e:   # the suite's reader refuses the frame; the library judges it on load
            line["doc_error"] = str(e)
        path = os.path.join(self.spool, "inputs", safe(sid) + ".mgt")
        with open(path, "wb") as fh:
            fh.write(data)
        try:
            msg = self.ask({"cmd": "call", "tool": "graphlet_load", "args": {"path": path}}, self.r.op_timeout)
        except Supervised as e:
            line["error"] = {"code": e.code, "fields": dict(e.fields, stage="load")}
            return
        answer = msg["result"].get("answer") if msg["ok"] else None
        if not msg["ok"]:
            line["error"] = {"code": "document_refused", "fields": {"stage": "load", "type": msg["type"]}}
            line["error_text"] = msg["error"]
        elif not isinstance(answer, dict) or isinstance(answer.get("error"), str):
            line["error"] = {"code": "document_refused", "fields": {"stage": "load", "code": (answer or {}).get(
                "error")}}
            line["error_text"] = (answer or {}).get("message")
        else:
            self.handles[sid] = answer["handle"]
            self.owner[answer["handle"]] = (sid, None)
            line["load"] = self.symbolise(answer, sid)

    # -- ops --------------------------------------------------------------------------------------------------
    def run_scenario(self, scenario):
        sid = scenario["id"]
        g = scenario["graphlet"]
        default = (g.get("ways") or {}).get("default")
        answers = {}
        for op in g.get("ops") or []:
            ways = op.get("ways", default)
            if ways is not None and PLAN_WAY[self.name] not in ways:
                answers[op["id"]] = None
                continue
            rec, view = self.run_op(sid, op, answers)
            answers[op["id"]] = view
            self.write_line(self.ops_file, rec)
            self.counts["ops"] += 1
            self.counts["answered" if rec["answer"] is not None else
                        "selector_error" if rec["selector_error"] else "error"] += 1

    def run_op(self, sid, op, answers):
        """One op in this way -> (its record, the answer later selectors read: graphlet_paging.selector_view of the
        record's answer, or None)."""
        started = time.time()
        tool = op["tool"]
        doc = self.document(sid)
        rec = {"run_id": self.r.run_id, "scenario": sid, "way": self.name, "op": op["id"], "tool": tool,
               "args_resolved": None, "selector_error": None, "pages": 0, "paging_complete": None, "answer": None,
               "error": None, "service_fields": None,
               "doc": None if doc["doc"] is None else dict((k, doc["doc"][k]) for k in ("body_digest",
                                                                                         "envelope_digest")),
               "library": self.library, "elapsed_ms": 0.0, "started_utc": iso(started)}
        try:
            args = self.graphlets(graphlet_select.resolve(op.get("args") or {}, answers), sid)
        except graphlet_select.SelectorError as e:
            rec["selector_error"] = str(e)
            return rec, None
        key = HANDLE_ARG.get(tool)
        if key and key not in args:
            args[key] = "<handle>"
        rec["args_resolved"] = dict((k, "<cursor>" if k == "cursor" and isinstance(v, str) else v)
                                    for k, v in args.items())
        if key and args[key] == "<handle>" and sid not in self.handles:
            rec["error"] = doc["error"]   # no graphlet to call the tool on
            return rec, None
        send = self.unsymbolise(args, sid)
        paging = op.get("paging", "none")
        # the continuation the op's own arguments give (a plan's cursor, a sequence's start): graphlet_paging's
        # runaway rule compares the first page's continuation with it
        given = (("start", send.get("start")) if tool == gp.SEQUENCE_TOOL else ("cursor", send.get("cursor"))) \
            if (send.get("start") if tool == gp.SEQUENCE_TOOL else send.get("cursor")) is not None else None
        timeout = op.get("timeout_s") or self.r.op_timeout
        pages, nxt, error = [], None, None
        while True:   # graphlet_paging's loop: follow next_call() while more(), stop at an error page
            call = dict(send)
            if nxt is not None:
                call[nxt[0]] = nxt[1]
            try:
                msg = self.ask({"cmd": "call", "tool": tool, "args": call}, timeout)
            except Supervised as e:
                msg = {"ok": False, "supervised": e.code, "fields": e.fields}
            sc.write_gz_json(os.path.join(self.dir, "raw", "%s.%s.%d.json.gz" % (safe(sid), safe(op["id"]),
                                                                                rec["pages"])),
                             {"run_id": self.r.run_id, "scenario": sid, "way": self.name, "op": op["id"],
                              "tool": tool, "page": rec["pages"], "args": call, "reply": msg})
            self.counts["calls"] += 1
            rec["pages"] += 1
            if "supervised" in msg:   # the call in flight when the process hung or died (restarted since)
                error = {"code": msg["supervised"], "fields": dict(msg["fields"], page=rec["pages"] - 1)}
                break
            if "non_finite" in msg:   # the answer holds a NaN or an infinity: never written as non-standard JSON
                error = {"code": "non_finite", "fields": {"pointer": msg["non_finite"], "page": rec["pages"] - 1}}
                break
            if not msg["ok"]:
                error = {"code": "exception", "fields": {"type": msg["type"]}}
                break
            if "unknown_tool" in msg["result"]:
                error = {"code": "unknown_tool", "fields": {}}
                break
            rec["elapsed_ms"] += msg["result"]["call_ms"]
            answer = msg["result"]["answer"]
            if gp.error_page(answer):
                error = gp.error_record(answer)
                break
            pages.append(answer)
            nxt = gp.next_call(tool, answer)
            if not gp.more(tool, paging, pages, given):
                break
        rec["elapsed_ms"] = round(rec["elapsed_ms"], 3)
        if error is None:
            error = gp.loop_error(tool, paging, pages, given)
        if error is None:
            try:
                answer = gp.assemble(tool, pages, nxt is None)
            except gp.AssemblyError as e:
                error = {"code": "assembly_error", "fields": {"reason": str(e)}}
        if error is not None:
            rec["error"] = self.symbolise(error, sid)
            return rec, None
        if tool in MINTING and isinstance(answer, dict) and isinstance(answer.get("handle"), str) \
                and answer["handle"] not in self.owner:
            self.minted[(sid, op["id"])] = answer["handle"]
            self.owner[answer["handle"]] = (sid, op["id"])
        rec["answer"] = self.symbolise(answer, sid)
        rec["paging_complete"] = gp.paging_complete(tool, pages[-1], answer)
        if tool in FILE_TOOLS:
            rec["file"] = self.keep_file(sid, op["id"], answer)
        return rec, gp.selector_view(tool, rec["answer"], pages[-1])

    def keep_file(self, sid, op, answer):
        """FORMAT v2 B3: the file an export or save wrote -> {name, bytes, sha256, raw}. Its text is kept gzipped at
        `raw` (relative to the way's directory, under raw/) with the handles and the way's directories written as
        symbols, as everywhere in the record, so that W5 and the compare stage read the same content from every way
        (byte for byte the file wherever it names no handle: all but json what=walk and a view's mgt); `bytes` is the
        file's own size, `sha256` that of the kept text."""
        path = answer.get("path") if isinstance(answer, dict) else None
        if not isinstance(path, str) or not os.path.isfile(path):
            return {"missing": self.symbolise(path, sid) if isinstance(path, str) else None}
        with open(path, "rb") as fh:
            data = fh.read()
        name = self.symbolise(os.path.basename(path), sid)
        ext = os.path.splitext(name)[1]
        rel = "raw/%s.%s.file%s.gz" % (safe(sid), safe(op), safe(ext) if ext else "")
        text = self.symbolise(data.decode("utf-8", "surrogateescape"), sid).encode("utf-8", "surrogateescape")
        with gzip.open(os.path.join(self.dir, rel), "wb") as fh:
            fh.write(text)
        return {"name": name, "bytes": len(data), "sha256": hashlib.sha256(text).hexdigest(), "raw": rel}


# ================================================================================================ the runner
class Runner(object):
    def __init__(self, args):
        self.args = args
        self.limits = read_limits(args.limits)
        self.python = args.python or sys.executable
        self.op_timeout = args.op_timeout
        self.ways = args.ways
        self.hashseeds = dict(DEFAULT_HASHSEEDS, **args.hashseeds)
        if len(set(self.hashseeds[w] for w in self.ways)) < len(self.ways):
            raise SystemExit("run_graphlet_local: the ways need distinct hash seeds, not %s"
                             % dict((w, self.hashseeds[w]) for w in self.ways))
        self.roots, self.commits = {}, {}
        for w in self.ways:
            drift = w == "drift"
            root = args.drift_library_root if drift else args.library_root
            if not root:
                raise SystemExit("run_graphlet_local: way %s needs %s" % (w, "--drift-library-root" if drift
                                                                          else "--library-root"))
            library_path(root)
            self.roots[w] = root
            self.commits[w] = library_commit(root, args.drift_library_commit if drift else args.library_commit)
        if args.doc:
            self.offline = True
            self.docs = dict(parse_doc(d) for d in args.doc)
            plan = sc.read_json(args.plan)
            self.scenarios = plan_scenarios(plan, [parse_doc(d)[0] for d in args.doc])
            self.out = os.path.abspath(args.out)
            self.run_id = args.run_id or os.path.basename(self.out)
        else:
            self.offline = False
            self.run_dir = os.path.abspath(args.run_dir)
            manifest = sc.read_json(args.manifest)
            if not isinstance(manifest, dict) or "scenarios" not in manifest:
                raise SystemExit("run_graphlet_local: %s is not a suite manifest" % args.manifest)
            self.scenarios = [s for s in manifest["scenarios"] if isinstance(s, dict)
                              and isinstance(s.get("graphlet"), dict)]
            self.out = os.path.abspath(args.out or args.run_dir)
            self.run_id = args.run_id or os.path.basename(self.run_dir)
        problems = plan_problems(self.scenarios)
        if problems:
            raise SystemExit("run_graphlet_local: the plan is malformed:\n  " + "\n  ".join(problems[:20]))
        self.by_id = dict((s["id"], s) for s in self.scenarios)
        if args.scenarios:
            unknown = [s for s in args.scenarios if s not in self.by_id]
            if unknown:
                raise SystemExit("run_graphlet_local: no graphlet scenario %s" % ", ".join(unknown))
            self.selected = [s for s in self.scenarios if s["id"] in set(args.scenarios)]
        else:
            self.selected = list(self.scenarios)

    def has_document(self, sid):
        return sid in self.docs if self.offline else sid in self.by_id

    def source(self, way, sid):
        """{kind: response | mgt, path, rel} of scenario |sid|'s document in way |way|."""
        if self.offline:
            if sid not in self.docs:   # a plan scenario no --doc names: its ops record no_document
                return {"kind": "mgt", "path": "", "rel": None}
            path = os.path.abspath(self.docs[sid])
            kind = "response" if re.search(r"\.json(\.gz)?\Z", path) else "mgt"
            return {"kind": kind, "path": path, "rel": path}
        if way in ("W1", "W2"):
            arm = WAY_SOURCE[way]
            rel = "arms/%s/%s" % (arm, sc.raw_name(sid, arm, self.args.repeat))
        else:
            rel = self.args.mgt_template.format(id=safe(sid), repeat=self.args.repeat)
        kind = "response" if re.search(r"\.json(\.gz)?\Z", rel) else "mgt"
        return {"kind": kind, "path": os.path.join(self.run_dir, rel), "rel": rel}

    def source_note(self, way):
        if self.offline:
            return {"offline": dict((k, os.path.abspath(v)) for k, v in sorted(self.docs.items()))}
        return {"run_dir": self.run_dir, "arm": WAY_SOURCE[way], "repeat": self.args.repeat,
                "template": None if way in ("W1", "W2") else self.args.mgt_template}

    def run(self):
        """Every way in turn; a way that fails (its process does not start or restarts too often) is reported and
        the next way runs. -> {way: its counts, or {"failed": why}}."""
        summary = {}
        for name in self.ways:
            try:
                way = Way(self, name)
            except WayFailed as e:
                summary[name] = {"failed": str(e)}
                print("run_graphlet_local: %s" % e, file=sys.stderr)
                continue
            state = "failed"
            try:
                for scenario in self.selected:
                    way.run_scenario(scenario)
                state = "done"
            except WayFailed as e:
                summary[name] = {"failed": str(e)}
                print("run_graphlet_local: %s" % e, file=sys.stderr)
            finally:
                way.close(state)
            if state != "done":
                continue
            summary[name] = dict(way.counts, seconds=round(time.time() - way.started, 2), dir=way.dir,
                                 restarts=len(way.restarts))
            if not self.args.quiet:
                c = way.counts
                print("[%s] %d scenarios, %d ops: %d answered, %d errors, %d selector errors, %d calls, %d restarts, "
                      "%.1f s; %s" % (name, len(self.selected), c["ops"], c["answered"], c["error"],
                                       c["selector_error"], c["calls"], len(way.restarts), summary[name]["seconds"],
                                       way.ops_path))
        return summary


def parse_doc(value):
    """--doc [ID=]PATH -> (ID, PATH); ID defaults to the file name without .gz, .json, .mgt, .graphlet."""
    if "=" in value and not os.path.exists(value):
        sid, path = value.split("=", 1)
        return sid, path
    name = os.path.basename(value)
    for ext in (".gz", ".json", ".mgt", ".graphlet"):
        if name.endswith(ext):
            name = name[:-len(ext)]
    return name, value


def plan_scenarios(plan, doc_ids):
    """An offline plan: a manifest ({"scenarios": [...]}), a list of scenarios, one scenario {"id", "graphlet"}, or a
    bare {"ops": [...]} (the first --doc's)."""
    if isinstance(plan, dict) and isinstance(plan.get("scenarios"), list):
        scenarios = plan["scenarios"]
    elif isinstance(plan, list):
        scenarios = plan
    elif isinstance(plan, dict) and isinstance(plan.get("graphlet"), dict):
        scenarios = [plan]
    elif isinstance(plan, dict) and isinstance(plan.get("ops"), list):
        scenarios = [{"id": doc_ids[0], "graphlet": plan}]
    else:
        raise SystemExit("run_graphlet_local: the plan is no manifest, scenario list, scenario or {\"ops\": [...]}")
    return [s for s in scenarios if isinstance(s, dict) and isinstance(s.get("graphlet"), dict)]


def plan_problems(scenarios):
    """What makes a plan unrunnable in any way (a malformed plan stops the run before any way starts)."""
    out, ids = [], set()
    for s in scenarios:
        sid = s.get("id")
        if not isinstance(sid, str) or sid in ids:
            out.append("scenario id %r: a string, unique" % (sid,))
        ids.add(sid)
        g = s["graphlet"]
        ways = (g.get("ways") or {}).get("default")
        if ways is not None and not isinstance(ways, list):
            out.append("%s: graphlet.ways.default is a list" % sid)
        ops = g.get("ops")
        if not isinstance(ops, list):
            out.append("%s: graphlet.ops is a list" % sid)
            continue
        seen = set()
        for op in ops:
            oid = op.get("id") if isinstance(op, dict) else None
            where = "%s op %r" % (sid, oid)
            if not isinstance(oid, str) or oid in seen:
                out.append("%s: an op id is a string, unique in its scenario" % where)
            seen.add(oid)
            if not isinstance(op, dict):
                continue
            if not isinstance(op.get("tool"), str):
                out.append("%s: tool is a tool name" % where)
            if not isinstance(op.get("args", {}), dict):
                out.append("%s: args is an object" % where)
            problem = gp.paging_problem(op.get("paging", "none"))
            if problem:
                out.append("%s: %s" % (where, problem))
            if "ways" in op and not isinstance(op["ways"], list):
                out.append("%s: ways is a list" % where)
            limit = op.get("timeout_s")
            if "timeout_s" in op and (isinstance(limit, bool) or not isinstance(limit, (int, float)) or limit <= 0):
                out.append("%s: timeout_s is a number of seconds > 0, not %r" % (where, limit))
    return out


# ================================================================================================ drift, export
DRIFT_FIELDS = ("args_resolved", "selector_error", "error", "answer", "paging_complete")


def drift_report(out, quiet=False):
    """The drift lane's report (DESIGN-v8 §2, owner D2): per (scenario, op) the record fields in which the branch
    head's library (arms/drift-lib) answers otherwise than the deployed commit's (arms/mcp-lib) on the same documents,
    an export's file by its digest. Report-only: written to arms/drift-lib/drift_report.json."""
    def records(way):
        rows, bad = sc.read_jsonl(os.path.join(out, "arms", WAY_DIRS[way], "graphlet_ops.jsonl"))
        return dict(((r["scenario"], r["op"]), r) for r in rows), bad
    deployed, bad_a = records("W4")
    head, bad_b = records("drift")
    if not deployed or not head:
        raise SystemExit("run_graphlet_local: --drift-report needs records of both W4 and drift under %s" % out)
    differ = []
    for key in sorted(set(deployed) & set(head)):
        a, b = deployed[key], head[key]
        what = [f for f in DRIFT_FIELDS if a.get(f) != b.get(f)]
        if (a.get("file") or {}).get("sha256") != (b.get("file") or {}).get("sha256"):
            what.append("file")
        if what:
            differ.append({"scenario": key[0], "op": key[1], "tool": a.get("tool"), "fields": what})
    report = {"deployed": sorted(set((r.get("library") or {}).get("commit") for r in deployed.values()), key=str),
              "head": sorted(set((r.get("library") or {}).get("commit") for r in head.values()), key=str),
              "compared": len(set(deployed) & set(head)), "differ": differ,
              "only_deployed": [list(k) for k in sorted(set(deployed) - set(head))],
              "only_head": [list(k) for k in sorted(set(head) - set(deployed))], "torn_lines": bad_a + bad_b,
              "written_utc": iso(time.time())}
    path = os.path.join(out, "arms", WAY_DIRS["drift"], "drift_report.json")
    sc.write_json_atomic(path, report)
    if not quiet:
        print("drift: %d ops compared, %d differ (%s vs %s); %s" % (report["compared"], len(differ),
                                                                    ",".join(map(str, report["deployed"])),
                                                                    ",".join(map(str, report["head"])), path))
    return report


def export_library(commit, to, repo=None, quiet=False):
    """`git archive COMMIT metagraph/api/python | tar -x -C TO`, and TO/LIBRARY_COMMIT naming the commit: read-only for
    the repository (no worktree, no checkout). A TO that already holds that commit is left as it is."""
    repo = repo or git(HERE, "rev-parse", "--show-toplevel")
    full = git(repo, "rev-parse", "--verify", commit + "^{commit}")
    marker = os.path.join(to, "LIBRARY_COMMIT")
    if os.path.isfile(marker):
        with open(marker, "r", encoding="utf-8") as fh:
            held = fh.read().strip()
        if held != full:
            raise SystemExit("run_graphlet_local: %s holds the library of %s, not %s" % (to, held, full))
        if not quiet:
            print("library %s already in %s" % (full, to))
        return full
    if os.path.exists(os.path.join(to, "metagraph")):
        raise SystemExit("run_graphlet_local: %s holds a metagraph/ without LIBRARY_COMMIT: use an empty directory"
                         % to)
    os.makedirs(to, exist_ok=True)
    archive = subprocess.Popen(["git", "-C", repo, "archive", full, "metagraph/api/python"], stdout=subprocess.PIPE)
    untar = subprocess.Popen(["tar", "-x", "-C", to], stdin=archive.stdout)
    archive.stdout.close()
    if untar.wait() != 0 or archive.wait() != 0:
        raise SystemExit("run_graphlet_local: git archive %s | tar -x -C %s failed" % (full, to))
    with open(marker, "w", encoding="utf-8") as fh:
        fh.write(full + "\n")
    if not quiet:
        print("library %s exported to %s" % (full, to))
    return full


# ================================================================================================ main
def parse_hashseeds(text):
    out = {}
    for item in (text or "").split(","):
        if not item:
            continue
        way, _, seed = item.partition("=")
        if way not in WAYS or not seed.isdigit() or not 1 <= int(seed) <= 4294967295:
            raise SystemExit("run_graphlet_local: --hashseed WAY=N, N in 1..4294967295 (0 turns hashing random "
                             "off), not %r" % item)
        out[way] = int(seed)
    return out


def parse_args(argv):
    ap = argparse.ArgumentParser(description="Run the traversal library over area G's documents and op plans "
                                             "(DESIGN-v8 §2; FORMAT-v8 §3).")
    ap.add_argument("--run-dir", help="the run directory holding the arms' raw answers (manifest mode)")
    ap.add_argument("--manifest", help="the suite manifest (its scenarios with a graphlet block)")
    ap.add_argument("--doc", action="append", help="offline: [ID=]PATH of a document (repeatable)")
    ap.add_argument("--plan", help="offline: the op plan (a manifest, scenarios, one scenario or {\"ops\": [...]})")
    ap.add_argument("--out", help="where arms/<way-dir>/ is written (default: --run-dir)")
    ap.add_argument("--ways", default=None, help="comma-separated ways of W1,W2,W4,drift (default W1,W2,W4; "
                                                  "offline W4)")
    ap.add_argument("--scenarios", default=None, help="comma-separated scenario ids (default: every graphlet one)")
    ap.add_argument("--library-root", help="a directory holding metagraph/api/python (W1, W2, W4)")
    ap.add_argument("--library-commit", help="the commit of --library-root when it records none")
    ap.add_argument("--drift-library-root", help="the branch head's library (the drift way)")
    ap.add_argument("--drift-library-commit")
    ap.add_argument("--limits", help="JSON {max_bytes, sequence_max_bytes, local_limits, max_graphlet_mb}")
    ap.add_argument("--repeat", type=int, default=0, help="the direct arms' repeat whose answers W1/W2 read (0)")
    ap.add_argument("--mgt-template", default=DEFAULT_MGT_TEMPLATE,
                    help="W4's and drift's document relative to --run-dir, with {id} (and {repeat}) (%s)"
                         % DEFAULT_MGT_TEMPLATE)
    ap.add_argument("--hashseed", default=None, help="WAY=N,... (defaults %s)" % ",".join(
        "%s=%d" % kv for kv in sorted(DEFAULT_HASHSEEDS.items())))
    ap.add_argument("--python", help="the interpreter of the ways (default: this one)")
    ap.add_argument("--op-timeout", type=float, default=DEFAULT_OP_TIMEOUT_S,
                    help="seconds a library call may take when its op states no timeout_s (%g); a call over it is "
                         "recorded `timeout` and the way's process restarted" % DEFAULT_OP_TIMEOUT_S)
    ap.add_argument("--run-id", default=None)
    ap.add_argument("--export-library", metavar="COMMIT", help="export COMMIT's metagraph/api/python to --to")
    ap.add_argument("--to", help="with --export-library: the (empty) target directory")
    ap.add_argument("--repo", help="with --export-library: the repository (default: this script's)")
    ap.add_argument("--drift-report", action="store_true", help="list the drift way's differences from W4 (--out)")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args(argv)
    if args.export_library or args.drift_report:
        if args.export_library and not args.to:
            ap.error("--export-library needs --to")
        if args.drift_report and not args.out:
            ap.error("--drift-report needs --out (the run directory)")
        return args
    if args.doc:
        if not args.plan or not args.out:
            ap.error("--doc needs --plan and --out")
    elif not args.run_dir or not args.manifest:
        ap.error("give --run-dir and --manifest, or --doc, --plan and --out")
    args.ways = [w for w in (args.ways or ("W4" if args.doc else "W1,W2,W4")).split(",") if w]
    bad = [w for w in args.ways if w not in WAYS]
    if bad or len(set(args.ways)) < len(args.ways):
        ap.error("--ways: distinct names of %s, not %s" % (", ".join(WAYS), ",".join(args.ways)))
    args.scenarios = [s for s in args.scenarios.split(",") if s] if args.scenarios else None
    args.hashseeds = parse_hashseeds(args.hashseed)
    return args


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if argv[:1] == [CHILD]:
        return child_main(argv)
    args = parse_args(argv)
    if args.export_library:
        export_library(args.export_library, args.to, args.repo, args.quiet)
        return 0
    if args.drift_report:
        drift_report(os.path.abspath(args.out), args.quiet)
        return 0
    summary = Runner(args).run()
    return 1 if any("failed" in v for v in summary.values()) else 0


if __name__ == "__main__":
    sys.exit(main())
