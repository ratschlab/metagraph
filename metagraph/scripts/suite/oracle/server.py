"""A `metagraph server_query` on a free 127.0.0.1 port, and a small stdlib HTTP client.

Only localhost is ever contacted. The server is stopped by stop() (SIGTERM, then SIGKILL after
10 s), which run_oracle.py calls in a finally block and from its signal handlers.
"""

import json
import os
import socket
import subprocess
import time
import urllib.error
import urllib.request

HOST = '127.0.0.1'


def free_port():
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind((HOST, 0))
        return s.getsockname()[1]


class ServerError(RuntimeError):
    pass


class Server:
    def __init__(self, metagraph, args, log_path, timeout=120):
        self.metagraph = metagraph
        self.args = list(args)
        self.log_path = log_path
        self.process = None
        self.port = None
        self._log = None
        for _ in range(10):
            self.port = free_port()
            self._log = open(log_path, 'ab')
            cmd = ([metagraph, 'server_query'] + self.args
                   + ['--port', str(self.port), '--address', HOST, '-p', '2'])
            self._log.write((' '.join(cmd) + '\n').encode())
            self._log.flush()
            self.process = subprocess.Popen(cmd, stdout=self._log, stderr=subprocess.STDOUT,
                                            stdin=subprocess.DEVNULL)
            if self._wait_ready(timeout):
                return
            self.stop()
        raise ServerError(f'could not start server_query (log: {log_path})')

    @property
    def pid(self):
        return self.process.pid if self.process else None

    def _wait_ready(self, timeout):
        deadline = time.time() + timeout
        while time.time() < deadline:
            if self.process.poll() is not None:
                return False
            status, _ = self.get('stats', timeout=5)
            if status == 200:
                return True
            time.sleep(0.2)
        return False

    def url(self, route):
        return f'http://{HOST}:{self.port}/{route}'

    def get(self, route, timeout=60):
        return _http('GET', self.url(route), None, timeout)

    def post(self, route, payload, timeout=300):
        return _http('POST', self.url(route), json.dumps(payload).encode(), timeout)

    def stop(self):
        if self.process is not None and self.process.poll() is None:
            self.process.terminate()
            try:
                self.process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                self.process.kill()
                self.process.wait()
        if self._log is not None and not self._log.closed:
            self._log.close()

    def alive(self):
        return self.process is not None and self.process.poll() is None


def _http(method, url, data, timeout):
    """-> (status, parsed JSON body or text); status 0 when no answer came."""
    req = urllib.request.Request(url, data=data, method=method,
                                 headers={'Content-Type': 'application/json'})
    # never through a proxy: the server is on localhost
    opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))
    try:
        with opener.open(req, timeout=timeout) as r:
            body = r.read()
            status = r.status
    except urllib.error.HTTPError as e:
        body = e.read()
        status = e.code
    except (urllib.error.URLError, OSError, socket.timeout):
        return 0, None
    try:
        return status, json.loads(body)
    except ValueError:
        return status, body.decode(errors='replace')
