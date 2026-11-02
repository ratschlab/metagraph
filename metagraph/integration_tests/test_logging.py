import os
import shlex
import subprocess
import tempfile
import unittest

from parameterized import parameterized

from base import METAGRAPH


class TestLogging(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tempdir.cleanup)
        self.fasta = os.path.join(self.tempdir.name, 'input.fa')
        with open(self.fasta, 'w') as out:
            out.write('>sequence\nACGTACGTACGTACGT\n')

    def run_build(self, level, verbose=None, k=3):
        env = os.environ.copy()
        env.pop('SPDLOG_LEVEL', None)
        if level is not None:
            env['SPDLOG_LEVEL'] = level
        command = shlex.split(METAGRAPH) + [
            'build', '-k', str(k), '--suffix', 'A',
            '-o', os.path.join(self.tempdir.name, 'graph'), self.fasta,
        ]
        if verbose:
            command.append(verbose)
        return subprocess.run(command, env=env, capture_output=True,
                              text=True, timeout=30)

    @parameterized.expand([
        ('unset', None, None, True, False),
        ('empty', '', None, True, False),
        ('info', 'info', None, True, False),
        ('warn', 'warn', None, False, False),
        ('uppercase_warn', 'WARN', None, False, False),
        ('off', 'off', None, False, False),
        ('trace', 'trace', None, True, True),
        ('invalid', 'invalid-level', None, True, False),
        ('verbose_default', None, '-v', True, True),
        ('verbose_warn', 'warn', '--verbose', True, True),
        ('verbose_off', 'off', '-v', True, True),
    ])
    def test_level_precedence(self, name, level, verbose, info, trace):
        result = self.run_build(level, verbose)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual('[info]' in result.stdout, info, result.stdout)
        self.assertEqual('[trace]' in result.stderr, trace, result.stderr)
        self.assertNotIn('[trace]', result.stdout)
        self.assertEqual('failed to parse SPDLOG_LEVEL' in result.stderr,
                         level == 'invalid-level', result.stderr)

    @parameterized.expand([
        ('warn', 'warn', True),
        ('error', 'error', True),
        ('off', 'off', False),
    ])
    def test_error_filtering_preserves_exit_status(self, name, level, error):
        result = self.run_build(level, k=1)
        self.assertEqual(result.returncode, 1, result.stderr)
        self.assertEqual('[error]' in result.stderr, error, result.stderr)
        self.assertNotIn('[info]', result.stdout)
        self.assertNotIn('[trace]', result.stderr)
