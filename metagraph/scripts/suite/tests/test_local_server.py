"""The localhost server_query clients of the harness (SD-026): oracle/server.py's, scripts/traversal's
pattern_fixtures.py and traverse_coords_snapshots.py: a request to the local server never goes through an HTTP
proxy the environment names, and the coordinate snapshots' launcher retries another port when its server exits."""

import importlib.util
import json
import os
import sys
import tempfile
import threading
import unittest
import urllib.request
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SUITE = os.path.dirname(HERE)
TRAVERSAL = os.path.normpath(os.path.join(SUITE, "..", "traversal"))
sys.path.insert(0, SUITE)

from oracle import server as oracle_server  # noqa: E402


def load(name):
    spec = importlib.util.spec_from_file_location(name + "_under_test", os.path.join(TRAVERSAL, name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class Handler(BaseHTTPRequestHandler):
    def _answer(self, body):
        data = json.dumps(body).encode()
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def do_GET(self):
        self._answer({"route": self.path})

    def do_POST(self):
        n = int(self.headers.get("Content-Length") or 0)
        self._answer({"route": self.path, "echo": json.loads(self.rfile.read(n) or b"null")})

    def log_message(self, *args):
        pass


# a stand-in `server_query` for the launcher: its first start exits at once (the port taken), later ones serve
FAKE = r'''import json, os, sys
from http.server import BaseHTTPRequestHandler, HTTPServer
marker = os.environ["FAKE_MARKER"]
first = not os.path.exists(marker)
open(marker, "a").write("x")
if first:
    sys.exit(1)
port = int(sys.argv[sys.argv.index("--port") + 1])

class H(BaseHTTPRequestHandler):
    def _ok(self, body):
        data = json.dumps(body).encode()
        self.send_response(200)
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def do_GET(self):
        self._ok({"route": self.path})

    def do_POST(self):
        n = int(self.headers.get("Content-Length") or 0)
        self.rfile.read(n)
        self._ok({"results": []})

    def log_message(self, *a):
        pass

HTTPServer(("127.0.0.1", port), H).serve_forever()
'''


class LocalServerClients(unittest.TestCase):
    def setUp(self):
        self.httpd = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
        self.port = self.httpd.server_address[1]
        threading.Thread(target=self.httpd.serve_forever, daemon=True).start()
        self.addCleanup(self.httpd.server_close)
        self.addCleanup(self.httpd.shutdown)
        # a proxy that answers nothing, and no NO_PROXY: a client that honours the environment's proxy fails
        env = dict((k, v) for k, v in os.environ.items() if k.lower() != "no_proxy")
        env.update(HTTP_PROXY="http://127.0.0.1:9", http_proxy="http://127.0.0.1:9")
        patcher = mock.patch.dict(os.environ, env, clear=True)
        patcher.start()
        self.addCleanup(patcher.stop)
        urllib.request._opener = None   # urlopen's default opener reads the proxies when it is built

    def test_no_client_goes_through_a_proxy(self):
        status, body = oracle_server._http("GET", "http://127.0.0.1:%d/stats" % self.port, None, 5)
        self.assertEqual((status, body), (200, {"route": "/stats"}))
        pf = load("pattern_fixtures")
        status, raw = pf.http(self.port, "GET", "/stats")
        self.assertEqual((status, json.loads(raw)), (200, {"route": "/stats"}))
        tc = load("traverse_coords_snapshots")
        s = object.__new__(tc._Server)
        s.url = "http://127.0.0.1:%d" % self.port
        status, raw = s.post({"seeds": []})
        self.assertEqual((status, json.loads(raw)["echo"]), (200, {"seeds": []}))

    def test_the_snapshots_launcher_retries_a_port(self):
        tc = load("traverse_coords_snapshots")
        with tempfile.TemporaryDirectory() as d:
            binary = os.path.join(d, "metagraph")
            with open(binary, "w") as fh:
                fh.write("#!%s\n%s" % (sys.executable, FAKE))
            os.chmod(binary, 0o755)
            os.environ["FAKE_MARKER"] = os.path.join(d, "starts")
            s = tc._Server(binary, d)
            try:
                status, raw = s.post({"seeds": []})
                self.assertEqual((status, json.loads(raw)), (200, {"results": []}))
            finally:
                s.close()
            with open(os.path.join(d, "starts")) as fh:
                self.assertEqual(fh.read(), "xx")


if __name__ == "__main__":
    unittest.main()
