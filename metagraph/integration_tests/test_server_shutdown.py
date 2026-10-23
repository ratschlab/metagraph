import signal
import socket
import subprocess
import tempfile
import time

import requests
from parameterized import parameterized

from base import METAGRAPH, TEST_DATA_DIR, TestingBase


class TestServerShutdown(TestingBase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.graph = cls.tempdir.name + '/graph.dbg'
        cls.annotation = cls.tempdir.name + '/annotation'
        fasta = TEST_DATA_DIR + '/transcripts_100.fa'
        cls._build_graph(fasta, cls.graph, 6, 'succinct')
        cls._annotate_graph(fasta, cls.graph, cls.annotation, 'column')

    @parameterized.expand([
        ('term_one_thread', signal.SIGTERM, 1, False),
        ('term_multiple_threads', signal.SIGTERM, 4, False),
        ('term_mmap', signal.SIGTERM, 2, True),
        ('interrupt', signal.SIGINT, 2, False),
    ])
    def test_shutdown(self, name, shutdown_signal, threads, mmap):
        with socket.socket() as reservation:
            reservation.bind(('127.0.0.1', 0))
            port = reservation.getsockname()[1]

        command = [METAGRAPH, 'server_query', '-i', self.graph,
                   '-a', self.annotation + '.column.annodbg',
                   '--address', '127.0.0.1', '--port', str(port), '-p', str(threads)]
        if mmap:
            command.append('--mmap')

        with tempfile.TemporaryFile(mode='w+') as log:
            process = subprocess.Popen(command, stdout=log, stderr=log)
            try:
                with requests.Session() as client:
                    client.trust_env = False
                    deadline = time.monotonic() + 30
                    while time.monotonic() < deadline:
                        self.assertIsNone(process.poll(), 'Server exited before readiness')
                        try:
                            response = client.get(f'http://127.0.0.1:{port}/stats', timeout=0.2)
                            if response.status_code == 200:
                                break
                        except requests.RequestException:
                            pass
                        time.sleep(0.05)
                    else:
                        self.fail('Server did not become ready')

                # An incomplete request must not keep the HTTP worker alive.
                with socket.create_connection(('127.0.0.1', port), timeout=2) as connection:
                    connection.sendall(b'GET /stats HTTP/1.1\r\nHost: localhost\r\n')
                    process.send_signal(shutdown_signal)
                    self.assertEqual(process.wait(timeout=10), 0)
                log.seek(0)
                self.assertIn('[Server] Shutting down...', log.read())
            finally:
                if process.poll() is None:
                    process.kill()
                process.wait(timeout=10)
