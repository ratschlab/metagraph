import errno
import os
import re
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
        command = [METAGRAPH, 'server_query', '-i', self.graph,
                   '-a', self.annotation + '.column.annodbg',
                   '--address', '127.0.0.1', '--port', '0', '-p', str(threads)]
        if mmap:
            command.append('--mmap')

        with tempfile.NamedTemporaryFile(mode='w+') as log:
            process = subprocess.Popen(command, stdout=log, stderr=log)
            try:
                with requests.Session() as client:
                    client.trust_env = False
                    deadline = time.monotonic() + 30
                    while time.monotonic() < deadline:
                        self.assertIsNone(process.poll(), 'Server exited before readiness')
                        with open(log.name) as startup_log:
                            match = re.search(r'\[Server\] Listening on port (\d+)',
                                              startup_log.read())
                        if not match:
                            time.sleep(0.05)
                            continue
                        port = int(match.group(1))
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

    @parameterized.expand([
        ('term', signal.SIGTERM),
        ('interrupt', signal.SIGINT),
    ])
    def test_shutdown_during_startup(self, name, shutdown_signal):
        # Hold startup at the index-list read, after signal registration but
        # before graph loading or HTTP startup. The FIFO synchronizes the signal.
        with tempfile.TemporaryDirectory() as tmp:
            index_list = os.path.join(tmp, 'indexes.csv')
            os.mkfifo(index_list)
            command = [METAGRAPH, 'server_query', index_list, '--address', '127.0.0.1',
                       '--port', '0', '-p', '2']
            with tempfile.TemporaryFile(mode='w+') as log:
                process = subprocess.Popen(command, stdout=log, stderr=log)
                writer = None
                try:
                    deadline = time.monotonic() + 30
                    while time.monotonic() < deadline:
                        self.assertIsNone(process.poll(), 'Server exited before opening index list')
                        try:
                            writer = os.open(index_list, os.O_WRONLY | os.O_NONBLOCK)
                            break
                        except OSError as error:
                            if error.errno != errno.ENXIO:
                                raise
                            time.sleep(0.05)
                    else:
                        self.fail('Server did not open the index list')

                    process.send_signal(shutdown_signal)
                    index = f'test,{self.graph},{self.annotation}.column.annodbg\n'.encode()
                    self.assertEqual(os.write(writer, index), len(index))
                    os.close(writer)
                    writer = None
                    self.assertEqual(process.wait(timeout=30), 0)
                    log.seek(0)
                    output = log.read()
                    self.assertIn('[Server] All graphs were loaded', output)
                    self.assertIn('[Server] Shutting down...', output)
                finally:
                    if writer is not None:
                        os.close(writer)
                    if process.poll() is None:
                        process.kill()
                    process.wait(timeout=10)
