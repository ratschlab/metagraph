"""A local stand-in for the staging host, for the runner's offline tests (stdlib http.server in a thread).

It serves the routes the runner touches under an optional path prefix (like the public URL's
/api/refseq33m-experimental): the three capabilities documents, POST /pattern (gzip when asked; 400 for an
unknown field), /resolve, /traverse, and /slow, which sleeps for the `sleep_s` of its body.
"""

import gzip
import json
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

RELEASE = "refseq97-b7214d29a934ff89512ad254bbb66ba793842928"
SERVER_INSTANCE = "feedfacecafebeef"


def capabilities():
    return {"release": RELEASE, "server_instance": SERVER_INSTANCE, "ready": True,
            "deadline_check": {"observed_max_uninterruptible_ms": int(time.time() * 1000) % 1000},
            "pattern": {"k": 31, "index_fp": None}}


def pattern_capabilities():
    return {"k": 31, "index": {"index_fp": "907ec986aaaa", "release": RELEASE}}


class Handler(BaseHTTPRequestHandler):
    prefix = ""
    log = None

    def log_message(self, fmt, *args):  # quiet
        pass

    def _path(self):
        path = self.path
        if self.prefix and path.startswith(self.prefix):
            path = path[len(self.prefix):]
        return path

    def _send(self, status, obj, compress=True):
        data = json.dumps(obj).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json")
        if compress and "gzip" in (self.headers.get("Accept-Encoding") or ""):
            data = gzip.compress(data)
            self.send_header("Content-Encoding", "gzip")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def do_GET(self):
        path = self._path()
        if self.log is not None:
            self.log.append(("GET", path))
        if path == "/capabilities":
            self._send(200, capabilities())
        elif path == "/pattern/capabilities":
            self._send(200, pattern_capabilities())
        elif path == "/traverse/capabilities":
            self._send(200, {"server_instance": SERVER_INSTANCE, "release": RELEASE})
        else:
            self._send(404, {"error": "not found"}, compress=False)

    def do_POST(self):
        path = self._path()
        length = int(self.headers.get("Content-Length") or 0)
        raw = self.rfile.read(length) if length else b""
        if self.log is not None:
            self.log.append(("POST", path))
        try:
            body = json.loads(raw.decode("utf-8")) if raw else None
        except ValueError:
            self._send(400, {"code": "invalid_request", "error": "request: not JSON"}, compress=False)
            return
        if path == "/pattern":
            if isinstance(body, dict) and "bogus" in body:
                self._send(400, {"code": "invalid_request", "error": "request: unknown field 'bogus'"},
                           compress=False)
                return
            pats = []
            for p in (body or {}).get("patterns", []):
                pats.append({"id": p.get("id"), "counts": {"contexts": {"value": 1, "relation": "exact"}},
                             "timing": {"elapsed_ms": 1.5},
                             "work": {"steps": 10, "ranges_visited": 10, "memory_bytes": 1000}})
            self._send(200, {"patterns": pats, "timing": {"elapsed_ms": 2.0}})
        elif path == "/resolve":
            self._send(200, {"k": 31, "num_kmers": 30, "labels": [], "timing": {"elapsed_ms": 3.0}})
        elif path == "/traverse":
            self._send(200, {"results": [{"arms": {"right": {"counters": {"steps": 5}},
                                                   "left": {"counters": {"steps": 7}}},
                                          "annotation": {"rows_requested": 4},
                                          "timing": {"walk_ms": 1.0}}]})
        elif path == "/slow":
            time.sleep(float((body or {}).get("sleep_s", 1.0)))
            self._send(200, {"slept": True})
        else:
            self._send(404, {"error": "not found"}, compress=False)


class FakeHost(object):
    def __init__(self, prefix="/api/refseq33m-experimental"):
        log = []

        class H(Handler):
            pass
        H.prefix = prefix
        H.log = log
        self.log = log
        self.server = ThreadingHTTPServer(("127.0.0.1", 0), H)
        self.server.daemon_threads = True
        self.port = self.server.server_address[1]
        self.base = "http://127.0.0.1:%d%s" % (self.port, prefix)
        self.thread = threading.Thread(target=self.server.serve_forever, daemon=True)

    def __enter__(self):
        self.thread.start()
        return self

    def __exit__(self, *exc):
        self.server.shutdown()
        self.server.server_close()
