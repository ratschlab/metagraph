"""frames(): pandas is imported lazily, by frames() alone; the package imports without
pandas (§5, §9)."""

import os
import subprocess
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import traverse_testlib as T  # noqa: E402

try:
    import pandas  # noqa: F401
    HAVE_PANDAS = True
except ImportError:
    HAVE_PANDAS = False


class TestImportWithoutPandas(unittest.TestCase):
    def test_package_imports_without_pandas(self):
        code = ("import sys; sys.modules['pandas'] = None; import metagraph.traverse as t; "
                "import metagraph.traverse.mcp_tools, metagraph.traverse.frames; "
                "g = t.parse(open(sys.argv[1]).read()); assert g.walks('right'); "
                "assert 'pandas' not in [m for m in sys.modules if sys.modules[m]]")
        p = subprocess.run([sys.executable, '-c', code, os.path.join(T.DOCS, 'merge.mgt')],
                           cwd=T.API, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(0, p.returncode, p.stderr.decode())

    def test_frames_needs_pandas_only_when_called(self):
        code = ("import sys; sys.modules['pandas'] = None; from metagraph.traverse import parse;"
                " from metagraph.traverse.frames import frames; "
                "g = parse(open(sys.argv[1]).read())\n"
                "try:\n    frames(g)\nexcept ImportError:\n    print('lazy')")
        p = subprocess.run([sys.executable, '-c', code, os.path.join(T.DOCS, 'merge.mgt')],
                           cwd=T.API, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(b'lazy\n', p.stdout, p.stderr.decode())


@unittest.skipUnless(HAVE_PANDAS, 'pandas is not installed')
class TestFrames(unittest.TestCase):
    def test_frames(self):
        from metagraph.traverse.frames import frames
        f = frames(T.body('merge'))
        self.assertEqual(7, len(f['segments']))
        # both.fa has its kept run and the two runs the merges closed; a.fa is displayed
        # from the merge at 62, b.fa from the one at 36 (each merge's first parent is the
        # one carried by the most labels)
        self.assertEqual(['c:0', 'c:1', 'c:2', 'c:3', 'c:3', 'c:3'], list(f['runs']['ref']))
        self.assertEqual([62, 36, 0, 0, 0, 0], list(f['runs']['evidence_from']))
        self.assertEqual(1, len(f['walks']))
        f = frames(T.body('annotate'))
        self.assertEqual(12, len(f['presence']))
        self.assertEqual(0, len(f['runs']))


if __name__ == '__main__':
    unittest.main()
