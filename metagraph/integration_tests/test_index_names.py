"""`metagraph stats --check-names`: the names an index carries, as the serving routes need them
(the owner's decision of 2026-10-11): every label name valid UTF-8, and, for a coordinate
annotation with its record mapping (.seqs) beside it, every record header valid UTF-8 and
distinct within its column. One JSON object alone on stdout, exit 0 when names_ok, else 1. And
scripts/traversal/index_manifest.py --metagraph recording that object under the manifest key
"names", outside the identity: index_fp is the same with and without it."""
import json
import os
import subprocess
import sys
import unittest
from subprocess import PIPE
from tempfile import TemporaryDirectory

from base import TestingBase, METAGRAPH, NUM_THREADS

REPO = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
MANIFEST_PY = os.path.join(REPO, 'scripts', 'traversal', 'index_manifest.py')

# the mini index (scripts/traversal/build_mini_refseq.sh), as test_traverse.py finds it
MINI_DIR = os.environ.get('METAGRAPH_MINI_REFSEQ', os.path.join(os.getcwd(), 'mini_refseq'))
MINI_GRAPH = 'graph_k31.dbg'
MINI_ANNO = 'annotation.relaxed.relabeled.row_diff_brwt_coord.annodbg'
MINI_FILES = (MINI_GRAPH, MINI_GRAPH + '.anchors', MINI_GRAPH + '.rd_succ', MINI_ANNO,
              'annotation.relaxed.relabeled.seqs')
_MINI_PRESENT = all(os.path.isfile(os.path.join(MINI_DIR, f)) for f in MINI_FILES)
_MINI_GUARD = (_MINI_PRESENT or os.environ.get('METAGRAPH_REQUIRE_GUARDS', '') == '1')

NAMES_OK = {'labels': 9, 'labels_invalid_utf8': 0, 'headers_checked': True, 'headers': 42,
            'headers_invalid_utf8': 0, 'duplicate_record_headers': 0, 'first_invalid': [],
            'names_ok': True}


def check_names(annotation):
    """(exit code, the object, stderr) of `stats --check-names -a ANNOTATION`: stdout holds the
    JSON object and nothing else."""
    res = subprocess.run([METAGRAPH, 'stats', '--check-names', '-a', annotation],
                         stdout=PIPE, stderr=PIPE)
    return (res.returncode, json.loads(res.stdout.decode('utf-8')),
            res.stderr.decode('utf-8', 'replace'))


class TestIndexNames(TestingBase):
    K = 11
    SEQ_A = 'ACGTTGCAGGCTAGCTTACGGATCCAGTACGATGCCTTAGC'
    SEQ_B = 'TTGACCGTAGGCATCGATTCAAGGCTTACGACTGGATACCA'
    SEQ_C = 'GGATCCTTAAGCGCGTATACCAGTTGCAATGGCTCGATCAG'
    SEQ_D = 'CATGCAGTTACGGATTCCGAATGGCTAGCTTAGCCGATACG'

    def setUp(self):
        self.tempdir = TemporaryDirectory()
        self.d = self.tempdir.name

    def _fasta(self, records, name='records.fa'):
        path = os.path.join(self.d, name)
        with open(path, 'wb') as f:
            for header, seq in records:
                f.write(b'>' + header + b'\n' + seq.encode() + b'\n')
        return path

    def _graph(self, fasta):
        self._build_graph(fasta, os.path.join(self.d, 'graph'), self.K, 'succinct')
        return os.path.join(self.d, 'graph.dbg')

    def _run(self, command):
        res = subprocess.run([command], shell=True, stdout=PIPE, stderr=PIPE)
        self.assertEqual(0, res.returncode, res.stderr.decode('utf-8', 'replace'))
        return res

    @unittest.skipUnless(_MINI_GUARD, "the mini index is not built (scripts/traversal/"
                                      "build_mini_refseq.sh)")
    def test_the_mini_index_names_are_ok(self):
        self.assertTrue(_MINI_PRESENT, f'the mini index is not in {MINI_DIR}')
        code, names, err = check_names(os.path.join(MINI_DIR, MINI_ANNO))
        self.assertEqual(0, code, err)
        self.assertEqual(NAMES_OK, names)

    def test_a_label_that_is_not_utf8_fails_the_check(self):
        # the utf8_labels fixture server's index (scripts/traversal/pattern_fixtures.py): a
        # header that is not UTF-8 and the header its replaced spelling would be
        fasta = self._fasta([(b'label\xff', self.SEQ_A), (b'label\xef\xbf\xbd', self.SEQ_B)])
        graph = self._graph(fasta)
        self._run(f'{METAGRAPH} annotate -p {NUM_THREADS} -i {graph} --anno-header '
                  f'--allow-unchecked-headers -o {self.d}/anno {fasta}')
        code, names, _ = check_names(f'{self.d}/anno.column.annodbg')
        self.assertEqual(1, code)
        self.assertEqual({'labels': 2, 'labels_invalid_utf8': 1, 'headers_checked': False,
                          'headers': 0, 'headers_invalid_utf8': 0, 'duplicate_record_headers': 0,
                          'first_invalid': [{'kind': 'label', 'problem': 'invalid_utf8',
                                             'column': 0, 'record': None, 'hex': '6c6162656cff'}],
                          'names_ok': False}, names)

    def test_check_names_takes_one_annotation_and_nothing_else(self):
        """Without an annotation, with two, or with a graph beside it, --check-names is refused
        on stderr with exit 1 and nothing on stdout, which holds the report alone."""
        fasta = self._fasta([(b'a', self.SEQ_A), (b'b', self.SEQ_B)])
        graph = self._graph(fasta)
        self._run(f'{METAGRAPH} annotate -p {NUM_THREADS} -i {graph} --anno-header '
                  f'-o {self.d}/anno {fasta}')
        anno = f'{self.d}/anno.column.annodbg'
        for args in ('', f'-a {anno} -a {anno}', f'{graph} -a {anno}'):
            res = subprocess.run([f'{METAGRAPH} stats --check-names {args}'], shell=True,
                                 stdout=PIPE, stderr=PIPE)
            self.assertEqual(1, res.returncode, args)
            self.assertEqual(b'', res.stdout, args)
            self.assertIn('--check-names takes exactly one annotation',
                          res.stderr.decode('utf-8', 'replace'), args)
        code, names, _ = check_names(anno)
        self.assertEqual(0, code)
        self.assertEqual(2, names['labels'])
        self.assertTrue(names['names_ok'])

    def test_a_repeated_record_header_fails_the_check(self):
        fasta = self._fasta([(b'dup', self.SEQ_A), (b'other', self.SEQ_B), (b'dup', self.SEQ_C),
                             (b'bad\xff', self.SEQ_D)])
        graph = self._graph(fasta)
        anno = f'{self.d}/anno'
        # a coordinate annotation with its record mapping, built unchecked
        self._run(f'{METAGRAPH} annotate -p {NUM_THREADS} -i {graph} --anno-filename '
                  f'--coordinates -o {anno} {fasta}')
        self._run(f'{METAGRAPH} transform_anno -p {NUM_THREADS} --anno-type column_coord '
                  f'--coordinates -o {anno} {anno}.column.annodbg')
        self._run(f'{METAGRAPH} annotate -p {NUM_THREADS} -i {graph} --anno-filename '
                  f'--index-header-coords --allow-unchecked-headers -o {anno} {fasta}')
        code, names, _ = check_names(f'{anno}.column_coord.annodbg')
        self.assertEqual(1, code)
        self.assertEqual({'labels': 1, 'labels_invalid_utf8': 0, 'headers_checked': True,
                          'headers': 4, 'headers_invalid_utf8': 1, 'duplicate_record_headers': 1,
                          'first_invalid': [{'kind': 'header', 'problem': 'invalid_utf8',
                                             'column': 0, 'record': 3, 'hex': '626164ff'},
                                            {'kind': 'header', 'problem': 'duplicate',
                                             'column': 0, 'record': 2, 'hex': '647570',
                                             'duplicate_of': 0}],
                          'names_ok': False}, names)
        # without the mapping only the labels are checked
        os.remove(f'{anno}.seqs')
        code, names, _ = check_names(f'{anno}.column_coord.annodbg')
        self.assertEqual(0, code)
        self.assertEqual({'labels': 1, 'labels_invalid_utf8': 0, 'headers_checked': False,
                          'headers': 0, 'headers_invalid_utf8': 0, 'duplicate_record_headers': 0,
                          'first_invalid': [], 'names_ok': True}, names)

    @unittest.skipUnless(_MINI_GUARD, "the mini index is not built (scripts/traversal/"
                                      "build_mini_refseq.sh)")
    def test_the_manifest_records_the_names_of_the_mini_index(self):
        self.assertTrue(_MINI_PRESENT, f'the mini index is not in {MINI_DIR}')
        graph, anno = os.path.join(MINI_DIR, MINI_GRAPH), os.path.join(MINI_DIR, MINI_ANNO)
        plain, named = f'{self.d}/plain.json', f'{self.d}/named.json'
        res = subprocess.run([sys.executable, MANIFEST_PY, '-i', graph, '-a', anno, '-o', plain],
                             stdout=PIPE, stderr=PIPE)
        self.assertEqual(0, res.returncode, res.stderr.decode())
        self.assertIn('names not checked (no --metagraph PATH)', res.stderr.decode())
        res = subprocess.run([sys.executable, MANIFEST_PY, '-i', graph, '-a', anno, '-o', named,
                              '--metagraph', METAGRAPH], stdout=PIPE, stderr=PIPE)
        self.assertEqual(0, res.returncode, res.stderr.decode())
        self.assertIn('names_ok true', res.stdout.decode())
        with open(plain) as f:
            p = json.load(f)
        with open(named) as f:
            n = json.load(f)
        self.assertNotIn('names', p)
        self.assertEqual(NAMES_OK, n['names'])
        self.assertEqual(p['index_fp'], n['index_fp'])
        self.assertEqual(p['files'], n['files'])
        # --verify re-checks against the bundle (the manifest lives apart from it: --root)
        res = subprocess.run([sys.executable, MANIFEST_PY, '--verify', named, '--root', MINI_DIR,
                              '--metagraph', METAGRAPH], stdout=PIPE, stderr=PIPE)
        self.assertEqual(0, res.returncode, res.stdout.decode() + res.stderr.decode())
        self.assertIn('names now: labels 9 (not UTF-8: 0), record headers 42 (not UTF-8: 0, '
                      'repeats: 0), names_ok true', res.stdout.decode())


if __name__ == '__main__':
    unittest.main()
