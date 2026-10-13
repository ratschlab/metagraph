"""scripts/traversal/index_manifest.py, batch mode: one manifest per (graph, annotation)
pair of a multi-graph server's list, each distinct file hashed once by parallel streams or
taken from precomputed sha256 digests, equal to what single mode writes; and the names check
(--metagraph, --check-names: `metagraph stats --check-names` recorded under "names", outside
the identity), driven here by a stand-in binary."""
import argparse
import contextlib
import hashlib
import io
import importlib.util
import json
import os
import sys
import tempfile
import threading
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, '..', '..', '..', 'scripts', 'traversal', 'index_manifest.py')


def _load():
    spec = importlib.util.spec_from_file_location('index_manifest', SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


im = _load()


def _args(**kw):
    # every option main() parses, at its default (--no-coord-mapping and --inventory too)
    a = dict(server_csv=None, jobs=1, digests=[], digests_only=False, match_base_names=False,
             out_dir=None, write_csv=None, force=False, graph=None, annotation=None, output=None,
             name=None, extra=[], no_coord_mapping=False, inventory=False, metagraph=None,
             check_names=None)
    a.update(kw)
    return argparse.Namespace(**a)


# what `metagraph stats --check-names` prints for an index whose names are all in order
NAMES_OK = {'labels': 3, 'labels_invalid_utf8': 0, 'headers': 5, 'headers_invalid_utf8': 0,
            'duplicate_record_headers': 0, 'first_invalid': [], 'names_ok': True}


class TestIndexManifestBatch(unittest.TestCase):
    def setUp(self):
        # the script reports to stdout and stderr; the tests read what it writes
        self.out, self.err = io.StringIO(), io.StringIO()
        self._quiet = contextlib.ExitStack()
        self._quiet.enter_context(contextlib.redirect_stdout(self.out))
        self._quiet.enter_context(contextlib.redirect_stderr(self.err))
        self.addCleanup(self._quiet.close)
        self.tmp = tempfile.TemporaryDirectory()
        d = self.root = self.tmp.name
        # two bundles sharing one graph (one graph, two annotations), and a third graph. The
        # coordinate annotation's .seqs is loaded for its own pair; the shared graph's mask and
        # Bloom filter are loaded for both of its pairs but are derived data, not part of any
        # identity; the files the server never opens for these pairs (row-diff anchors without
        # a .row_diff annotation, weights, a column's .coords) are in no bundle either: none of
        # these is ever hashed
        self.files = {}
        self.unloaded = {'g1.dbg.anchors', 'g1.dbg.weights', 'a2.column.annodbg.coords',
                         'g1.edgemask', 'g1.bloom'}
        for name, data in (('g1.dbg', b'graph one' * 1000), ('g1.edgemask', b'mask'),
                           ('g1.bloom', b'bloom'), ('g1.dbg.anchors', b'anchors'),
                           ('g1.dbg.weights', b'weights'),
                           ('a1.relaxed.row_diff_brwt_coord.annodbg', b'annotation one'),
                           ('a1.relaxed.seqs', b'headers'),
                           ('a2.column.annodbg', b'annotation two'),
                           ('a2.column.annodbg.coords', b'coords'),
                           ('g3.dbg', b'graph three'), ('a3.column.annodbg', b'annotation 3')):
            p = os.path.join(d, name)
            with open(p, 'wb') as f:
                f.write(data)
            self.files[name] = p
        f = self.files
        self.csv = os.path.join(d, 'graphs.csv')
        with open(self.csv, 'w') as out:
            out.write('\n'.join([
                'A,%s,%s,,ns.a' % (f['g1.dbg'], f['a1.relaxed.row_diff_brwt_coord.annodbg']),
                'B,%s,%s' % (f['g1.dbg'], f['a2.column.annodbg']),
                'B2,%s,%s' % (f['g1.dbg'], f['a2.column.annodbg']),   # the same pair again
                '',
                'C,%s,%s,%s' % (f['g3.dbg'], f['a3.column.annodbg'],
                                os.path.join(d, 'col4', 'c.manifest.json')),
            ]) + '\n')

    def tearDown(self):
        self.tmp.cleanup()

    def _single(self, graph, annotation, out, **kw):
        im.write(_args(graph=graph, annotation=annotation, output=out, **kw))
        with open(out) as f:
            return json.load(f)

    def _fake_metagraph(self):
        """A stand-in binary: `stats --check-names -a ANNOTATION` prints <ANNOTATION>.names.json
        and exits 1 when its names_ok is false, as the real one does; anything else exits 2."""
        path = os.path.join(self.root, 'metagraph')
        with open(path, 'w') as f:
            f.write('#!/bin/sh\n'
                    'test "$1" = stats -a "$2" = --check-names -a "$3" = -a || exit 2\n'
                    'cat "$4.names.json" || exit 2\n'
                    'grep -q \'"names_ok": true\' "$4.names.json"\n')
        os.chmod(path, 0o755)
        return path

    def _names(self, annotation, **kw):
        names = dict(NAMES_OK, **kw)
        with open(annotation + '.names.json', 'w') as f:
            json.dump(names, f, indent=1)
        return names

    def test_names_are_recorded_outside_the_identity(self):
        """--metagraph: the check's JSON is recorded under "names"; index_fp is the digest of
        the file list, the same with and without it; a failed check is recorded, and warned
        about, not refused."""
        f = self.files
        g, a = f['g1.dbg'], f['a1.relaxed.row_diff_brwt_coord.annodbg']
        binary = self._fake_metagraph()
        names = self._names(a)
        plain = self._single(g, a, os.path.join(self.root, 'plain.json'))
        self.assertNotIn('names', plain)
        self.assertIn('names not checked (no --metagraph PATH)', self.err.getvalue())
        got = self._single(g, a, os.path.join(self.root, 'named.json'), metagraph=binary)
        self.assertEqual(names, got['names'])
        self.assertEqual(plain['index_fp'], got['index_fp'])
        self.assertEqual(plain['files'], got['files'])
        self.assertIn('names_ok true', self.out.getvalue())
        self.assertNotIn('WARNING', self.err.getvalue())
        bad = self._names(a, labels_invalid_utf8=1, names_ok=False,
                          first_invalid=[{'kind': 'label', 'column': 2, 'record': None,
                                          'hex': '6c6162656cff', 'problem': 'invalid_utf8'}])
        got = self._single(g, a, os.path.join(self.root, 'bad.json'), metagraph=binary)
        self.assertEqual(bad, got['names'])
        self.assertEqual(plain['index_fp'], got['index_fp'])
        self.assertIn('WARNING', self.err.getvalue())
        self.assertIn('names_ok false', self.err.getvalue())

    def test_check_names_needs_the_binary_and_is_skipped_on_request(self):
        f = self.files
        g, a = f['g1.dbg'], f['a1.relaxed.row_diff_brwt_coord.annodbg']
        out = os.path.join(self.root, 'm.json')
        with self.assertRaises(im.ManifestError) as cm:
            im.write(_args(graph=g, annotation=a, output=out, check_names=True))
        self.assertIn('--metagraph', str(cm.exception))
        self.assertFalse(os.path.exists(out))
        binary = self._fake_metagraph()
        self._names(a)
        got = self._single(g, a, out, metagraph=binary, check_names=False)
        self.assertNotIn('names', got)
        # a binary that fails for another reason, or prints something else, is an error
        os.remove(a + '.names.json')
        with self.assertRaises(im.ManifestError):
            im.write(_args(graph=g, annotation=a, output=out, metagraph=binary))
        with open(a + '.names.json', 'w') as fh:
            fh.write('{"labels": 1, "names_ok": true}\n')
        with self.assertRaises(im.ManifestError) as cm:
            im.write(_args(graph=g, annotation=a, output=out, metagraph=binary))
        self.assertIn('headers', str(cm.exception))
        with self.assertRaises(im.ManifestError):
            im.write(_args(graph=g, annotation=a, output=out, metagraph=self.root + '/none'))

    def test_batch_records_names_per_pair(self):
        f = self.files
        binary = self._fake_metagraph()
        expected = {
            f['a1.relaxed.row_diff_brwt_coord.annodbg']:
                self._names(f['a1.relaxed.row_diff_brwt_coord.annodbg'], labels=1),
            f['a2.column.annodbg']: self._names(f['a2.column.annodbg'], labels=2, headers=0),
            f['a3.column.annodbg']: self._names(f['a3.column.annodbg'], labels=3,
                                                duplicate_record_headers=2, names_ok=False),
        }
        written = im.batch(_args(server_csv=self.csv, out_dir=self.root + '/b',
                                 metagraph=binary))
        self.assertEqual(3, len(written))
        for (g, a), out in written.items():
            with open(out) as fh:
                got = json.load(fh)
            self.assertEqual(expected[a], got['names'], out)
            ref = self._single(g, a, os.path.join(self.root, 'ref.json'))
            self.assertEqual(ref['index_fp'], got['index_fp'], out)
        self.assertIn('WARNING', self.err.getvalue())
        # skipped: no names (--force: line 5 fixes its pair's manifest path)
        written = im.batch(_args(server_csv=self.csv, out_dir=self.root + '/n', force=True))
        for out in written.values():
            with open(out) as fh:
                self.assertNotIn('names', json.load(fh))

    def test_verify_reports_and_compares_names(self):
        f = self.files
        g, a = f['g1.dbg'], f['a1.relaxed.row_diff_brwt_coord.annodbg']
        binary = self._fake_metagraph()
        self._names(a)
        out = os.path.join(self.root, 'v.json')
        self._single(g, a, out, metagraph=binary)
        im.verify(out)
        self.assertIn('names: labels 3 (not UTF-8: 0), record headers 5 (not UTF-8: 0, '
                      'repeats: 0), names_ok true', self.out.getvalue())
        self.assertNotIn('names now', self.out.getvalue())
        im.verify(out, metagraph=binary)
        self.assertIn('names now: labels 3', self.out.getvalue())
        # the names changed since the manifest was written: a problem
        self._names(a, headers_invalid_utf8=1, names_ok=False)
        with self.assertRaises(SystemExit):
            im.verify(out, metagraph=binary)
        self.assertIn('names differ from the recorded ones', self.out.getvalue())
        # a manifest without names: said so; the fresh check is reported, not a difference
        plain = os.path.join(self.root, 'p.json')
        self._single(g, a, plain)
        im.verify(plain, metagraph=binary)
        self.assertIn('names: not recorded', self.out.getvalue())
        self.assertIn('names now: labels 3 (not UTF-8: 0), record headers 5 (not UTF-8: 1',
                      self.out.getvalue())

    def _never(self, path, progress=None):
        raise AssertionError('%s was hashed before the names check was settled' % path)

    def test_names_refusals_come_before_any_hashing(self):
        """--check-names without a binary, a binary that cannot run and the binary's own
        failure are refused before a single file is hashed, in write and in batch mode (a
        refusal after hours of hashing would cost a full rerun), with nothing written."""
        f = self.files
        g, a = f['g1.dbg'], f['a1.relaxed.row_diff_brwt_coord.annodbg']
        out = os.path.join(self.root, 'm.json')
        not_executable = os.path.join(self.root, 'plain.txt')
        with open(not_executable, 'w') as fh:
            fh.write('not a binary\n')
        binary = self._fake_metagraph()     # no names.json beside the annotations: exits 2
        refusals = ((dict(check_names=True), '--metagraph'),
                    (dict(metagraph=self.root + '/none'), 'cannot run'),
                    (dict(metagraph=not_executable), 'cannot run'),
                    (dict(metagraph=binary), 'exited 2'))
        for kw, fragment in refusals:
            with self.assertRaises(im.ManifestError) as cm:
                im.write(_args(graph=g, annotation=a, output=out, **kw), hasher=self._never)
            self.assertIn(fragment, str(cm.exception))
            self.assertFalse(os.path.exists(out))
        csv_out = os.path.join(self.root, 'with.csv')
        for kw, fragment in refusals:
            with self.assertRaises(im.ManifestError) as cm:
                im.batch(_args(server_csv=self.csv, out_dir=self.root + '/b', write_csv=csv_out,
                               **kw), hasher=self._never)
            self.assertIn(fragment, str(cm.exception))
        self.assertFalse(os.path.exists(self.root + '/b'))
        self.assertFalse(os.path.exists(csv_out))

    def test_batch_refuses_one_pairs_failing_binary_before_hashing(self):
        """The binary failing on one pair (here: exit 2) refuses the whole batch before the
        hashing pool, naming the pair's line: no manifest of another pair is written and no
        digest is lost, as with every other refusal of the plan."""
        f = self.files
        binary = self._fake_metagraph()
        self._names(f['a1.relaxed.row_diff_brwt_coord.annodbg'])
        self._names(f['a3.column.annodbg'])     # a2's pair (line 2) has no names.json
        csv_out = os.path.join(self.root, 'with.csv')
        with self.assertRaises(im.ManifestError) as cm:
            im.batch(_args(server_csv=self.csv, out_dir=self.root + '/b', write_csv=csv_out,
                           metagraph=binary), hasher=self._never)
        self.assertIn('line 2', str(cm.exception))
        self.assertIn('exited 2', str(cm.exception))
        self.assertFalse(os.path.exists(self.root + '/b'))
        self.assertFalse(os.path.exists(os.path.join(self.root, 'col4')))
        self.assertFalse(os.path.exists(csv_out))
        # with every pair's names in order the batch is written, each pair's names recorded
        self._names(f['a2.column.annodbg'], labels=2)
        written = im.batch(_args(server_csv=self.csv, out_dir=self.root + '/b', write_csv=csv_out,
                                 metagraph=binary))
        self.assertEqual(3, len(written))
        for (g, a), out in written.items():
            with open(out) as fh:
                self.assertEqual(2 if a == f['a2.column.annodbg'] else 3,
                                 json.load(fh)['names']['labels'])
        self.assertTrue(os.path.exists(csv_out))

    def test_the_binarys_own_failure_reports_its_stderr(self):
        """`stats --check-names` exits 1 with nothing on stdout and the cause on stderr when it
        cannot read the annotation or its .seqs: that cause is reported, not "not JSON"."""
        f = self.files
        g, a = f['g1.dbg'], f['a1.relaxed.row_diff_brwt_coord.annodbg']
        out = os.path.join(self.root, 'm.json')
        failing = os.path.join(self.root, 'failing')
        with open(failing, 'w') as fh:
            fh.write('#!/bin/sh\necho "[error] Cannot read the labels of $4: truncated" >&2\n'
                     'exit 1\n')
        os.chmod(failing, 0o755)
        with self.assertRaises(im.ManifestError) as cm:
            im.write(_args(graph=g, annotation=a, output=out, metagraph=failing))
        self.assertIn('exited 1', str(cm.exception))
        self.assertIn('Cannot read the labels of ' + a + ': truncated', str(cm.exception))
        self.assertNotIn('not JSON', str(cm.exception))
        # output that is not the object: its stderr is reported with the parse error
        garbled = os.path.join(self.root, 'garbled')
        with open(garbled, 'w') as fh:
            fh.write('#!/bin/sh\necho "[info] loading"\necho "[warning] odd" >&2\nexit 0\n')
        os.chmod(garbled, 0o755)
        with self.assertRaises(im.ManifestError) as cm:
            im.write(_args(graph=g, annotation=a, output=out, metagraph=garbled))
        self.assertIn('not JSON', str(cm.exception))
        self.assertIn('[warning] odd', str(cm.exception))

    def test_the_skip_note_names_its_reason(self):
        f = self.files
        g, a = f['g1.dbg'], f['a1.relaxed.row_diff_brwt_coord.annodbg']
        binary = self._fake_metagraph()
        self._names(a)
        self._single(g, a, os.path.join(self.root, 'p.json'))
        self.assertIn('names not checked (no --metagraph PATH)', self.err.getvalue())
        self.assertNotIn('--no-check-names', self.err.getvalue())
        self.err.truncate(0)
        self.err.seek(0)
        got = self._single(g, a, os.path.join(self.root, 'n.json'), metagraph=binary,
                           check_names=False)
        self.assertNotIn('names', got)
        self.assertIn('names not checked (--no-check-names)', self.err.getvalue())
        self.assertNotIn('no --metagraph PATH', self.err.getvalue())
        self.err.truncate(0)
        self.err.seek(0)
        im.batch(_args(server_csv=self.csv, out_dir=self.root + '/b', metagraph=binary,
                       check_names=False))
        self.assertIn('names not checked (--no-check-names)', self.err.getvalue())
        self.assertNotIn('no --metagraph PATH', self.err.getvalue())

    def test_verify_honours_the_check_names_flags(self):
        """main(): --verify with --metagraph re-checks unless --no-check-names; --check-names
        without a binary, and a binary that cannot run, are refused in verify mode too."""
        f = self.files
        g, a = f['g1.dbg'], f['a1.relaxed.row_diff_brwt_coord.annodbg']
        binary = self._fake_metagraph()
        self._names(a)
        out = os.path.join(self.root, 'v.json')
        self._single(g, a, out, metagraph=binary)
        # the names changed since the manifest was written
        self._names(a, headers_invalid_utf8=1, names_ok=False)

        def run_main(*argv):
            with mock.patch.object(sys, 'argv', ['index_manifest.py'] + list(argv)):
                try:
                    im.main()
                except SystemExit as e:
                    return e.code
            return None

        self.assertIsNone(run_main('--verify', out, '--metagraph', binary, '--no-check-names'))
        self.assertIn('OK, index_fp', self.out.getvalue())
        self.assertNotIn('names now', self.out.getvalue())
        self.assertEqual(1, run_main('--verify', out, '--metagraph', binary))
        self.assertIn('names differ from the recorded ones', self.out.getvalue())
        code = run_main('--verify', out, '--check-names')
        self.assertIn('--metagraph', str(code))
        code = run_main('--verify', out, '--metagraph', self.root + '/none')
        self.assertIn('cannot run', str(code))

    def test_inventory_runs_no_names_check_and_probes_no_binary(self):
        """main(): --inventory prints the loader inventory and runs no check, so the names
        arguments are not probed there (a binary that cannot run, --check-names without one),
        while a write with the same arguments is refused before any file is hashed."""
        f = self.files
        g, a = f['g1.dbg'], f['a2.column.annodbg']

        def run_main(*argv):
            with mock.patch.object(sys, 'argv', ['index_manifest.py'] + list(argv)):
                try:
                    im.main()
                except SystemExit as e:
                    return e.code
            return None

        for extra in (('--metagraph', self.root + '/none'), ('--check-names',)):
            self.out.seek(0)
            self.out.truncate(0)
            self.assertIsNone(run_main('-i', g, '-a', a, '--inventory', *extra), extra)
            inventory = json.loads(self.out.getvalue())
            self.assertEqual([g, a], [e['path'] for e in inventory['files']], extra)
            self.assertIn('derived', inventory, extra)
        out = os.path.join(self.root, 'inv.json')
        code = run_main('-i', g, '-a', a, '-o', out, '--metagraph', self.root + '/none')
        self.assertIn('cannot run', str(code))
        self.assertFalse(os.path.exists(out))

    def test_batch_equals_single_mode_and_hashes_each_file_once(self):
        calls = []
        lock = threading.Lock()

        def hasher(path):
            with lock:
                calls.append(os.path.realpath(path))
            return im.sha256_file(path)

        out_dir = os.path.join(self.root, 'out')
        written = im.batch(_args(server_csv=self.csv, jobs=4, out_dir=out_dir,
                                 write_csv=os.path.join(self.root, 'with.csv')), hasher=hasher)
        self.assertEqual(3, len(written))
        # every distinct file once: the shared graph (and its sidecars) too; the files no
        # pair's server loads never
        self.assertEqual(sorted(set(calls)), sorted(calls))
        self.assertEqual(sorted(os.path.realpath(p) for n, p in self.files.items()
                                if n not in self.unloaded), sorted(calls))
        f = self.files
        for (g, a), out in written.items():
            with open(out) as fh:
                got = json.load(fh)
            ref = self._single(g, a, os.path.join(self.root, 'ref.json'))
            self.assertEqual(ref['index_fp'], got['index_fp'], out)
            self.assertEqual([(e['path'], e['size'], e['sha256']) for e in ref['files']],
                             [(e['path'], e['size'], e['sha256']) for e in got['files']])
            self.assertTrue(all(e['digest'] == 'hashed' for e in got['files']))
        # where they went: column 4, else the out dir by the first name
        self.assertEqual(os.path.join(self.root, 'col4', 'c.manifest.json'),
                         written[(f['g3.dbg'], f['a3.column.annodbg'])])
        self.assertEqual(os.path.join(out_dir, 'A.a1.relaxed.manifest.json'),
                         written[(f['g1.dbg'], f['a1.relaxed.row_diff_brwt_coord.annodbg'])])
        self.assertEqual(os.path.join(out_dir, 'B.a2.column.manifest.json'),
                         written[(f['g1.dbg'], f['a2.column.annodbg'])])
        with open(written[(f['g1.dbg'], f['a1.relaxed.row_diff_brwt_coord.annodbg'])]) as fh:
            self.assertEqual('ns.a', json.load(fh)['index_ns'])
        # the list written back: every line, the manifest column filled, index_ns kept, and
        # read again as the server reads it
        lines = im.parse_server_csv(os.path.join(self.root, 'with.csv'))
        self.assertEqual(['A', 'B', 'B2', 'C'], [c[0] for _, c in lines])
        for _, cols in lines:
            self.assertEqual(written[(cols[1], cols[2])], cols[3])
        self.assertEqual('ns.a', lines[0][1][4])
        self.assertEqual('', lines[1][1][4])
        # existing manifests are kept unless --force
        with self.assertRaises(im.ManifestError):
            im.batch(_args(server_csv=self.csv, out_dir=out_dir))
        im.batch(_args(server_csv=self.csv, out_dir=out_dir, force=True))

    def test_parallel_streams_write_the_same_manifests(self):
        a = im.batch(_args(server_csv=self.csv, jobs=1, out_dir=self.root + '/j1'))
        texts1 = {k: open(v).read() for k, v in a.items() if '/j1/' in v}
        b = im.batch(_args(server_csv=self.csv, jobs=4, out_dir=self.root + '/j4', force=True))
        texts4 = {k: open(v).read() for k, v in b.items() if '/j4/' in v}
        self.assertEqual(texts1, texts4)
        self.assertEqual(2, len(texts1))

    def test_precomputed_digests_replace_hashing(self):
        digests = os.path.join(self.root, 'sums.txt')
        with open(digests, 'w') as out:
            out.write('# from the storage checksums\n')
            for name, p in self.files.items():
                d = hashlib.sha256(open(p, 'rb').read()).hexdigest()
                # by exact path, by binary-mode marker, and by base name alone
                if name.startswith('g1'):
                    out.write('%s  %s\n' % (d, p))
                elif name.startswith('a1'):
                    out.write('%s *%s\n' % (d, p))
                else:
                    out.write('%s  /elsewhere/%s\n' % (d, name))

        def never(path):
            raise AssertionError('hashed %s' % path)

        # the lines by base name only match with --match-base-names
        with self.assertRaises(im.ManifestError):
            im.batch(_args(server_csv=self.csv, out_dir=self.root + '/o0', digests=[digests],
                           digests_only=True), hasher=never)
        written = im.batch(_args(server_csv=self.csv, out_dir=self.root + '/o', digests=[digests],
                                 digests_only=True, match_base_names=True), hasher=never)
        for (g, a), out in written.items():
            got = json.load(open(out))
            for e in got['files']:
                by_path = e['path'].startswith(('g1', 'a1'))
                self.assertEqual('precomputed' if by_path else 'precomputed-by-base-name',
                                 e['digest'], e['path'])
            ref = self._single(g, a, os.path.join(self.root, 'ref.json'))
            self.assertEqual(ref['index_fp'], got['index_fp'])

    def test_digest_problems_are_refused(self):
        p = self.files['g3.dbg']
        d = hashlib.sha256(open(p, 'rb').read()).hexdigest()
        # a file without a digest under --digests-only
        partial = os.path.join(self.root, 'partial.txt')
        with open(partial, 'w') as out:
            out.write('%s  %s\n' % (d, p))
        with self.assertRaises(im.ManifestError):
            im.batch(_args(server_csv=self.csv, out_dir=self.root + '/o', digests=[partial],
                           digests_only=True))
        # two digests for one path
        two = os.path.join(self.root, 'two.txt')
        with open(two, 'w') as out:
            out.write('%s  %s\n%s  %s\n' % (d, p, '0' * 64, p))
        with self.assertRaises(im.ManifestError):
            im.read_digests([two])
        # an ambiguous base name matches nothing: the file is hashed instead
        amb = os.path.join(self.root, 'amb.txt')
        with open(amb, 'w') as out:
            out.write('%s  /x/g3.dbg\n%s  /y/g3.dbg\n' % (d, '1' * 64))
        self.assertIsNone(im.DigestIndex(im.read_digests([amb]), base_names=True).find(p))
        # a malformed line
        bad = os.path.join(self.root, 'bad.txt')
        with open(bad, 'w') as out:
            out.write('xyz  %s\n' % p)
        with self.assertRaises(im.ManifestError):
            im.read_digests([bad])

    def test_digests_never_name_another_bundles_file(self):
        """Bundles share base names (graph.dbg). A checksum file written in one bundle's
        directory with relative paths must not name the other bundle's files by base name
        (their manifest would take the wrong digests), and a line spelling a path exactly must
        not skip the check for a second, different digest of the same file."""
        d = self.root
        for b, content in (('refseq', b'refseq graph'), ('uhgg', b'uhgg graph!!')):
            os.makedirs(os.path.join(d, b))
            for name, data in (('graph.dbg', content), ('annotation.column.annodbg', content * 2)):
                with open(os.path.join(d, b, name), 'wb') as f:
                    f.write(data)
        sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
        # sha256sum run inside refseq/, the file kept beside the bundles
        sums = os.path.join(d, 'refseq.sha256')
        with open(sums, 'w') as out:
            for name in ('graph.dbg', 'annotation.column.annodbg'):
                out.write('%s  %s\n' % (sha(os.path.join(d, 'refseq', name)), name))
        csv = os.path.join(d, 'two.csv')
        with open(csv, 'w') as out:
            for b in ('refseq', 'uhgg'):
                out.write('%s,%s/%s/graph.dbg,%s/%s/annotation.column.annodbg\n' % (b, d, b, d, b))
        # relative paths are the digest file's directory's: they name no bundle file here
        written = im.batch(_args(server_csv=csv, out_dir=d + '/o1'))
        for (g, a), out in written.items():
            got = json.load(open(out))
            self.assertTrue(all(e['digest'] == 'hashed' for e in got['files']))
            self.assertEqual(self._single(g, a, d + '/ref.json')['index_fp'], got['index_fp'])
        written = im.batch(_args(server_csv=csv, out_dir=d + '/o2', digests=[sums]))
        for (g, a), out in written.items():
            self.assertEqual(self._single(g, a, d + '/ref.json')['index_fp'],
                             json.load(open(out))['index_fp'])
        # by base name: one digest line for two different files is refused
        with self.assertRaises(im.ManifestError) as cm:
            im.batch(_args(server_csv=csv, out_dir=d + '/o3', digests=[sums],
                           match_base_names=True))
        self.assertIn('two different files', str(cm.exception))
        # the same file written inside refseq/: its lines name refseq's files by path, and are
        # never lent to uhgg's by base name, whatever the order of lookups
        inside = os.path.join(d, 'refseq', 'SHA256SUMS')
        with open(inside, 'w') as out:
            out.write(open(sums).read())
        written = im.batch(_args(server_csv=csv, out_dir=d + '/o4', digests=[inside],
                                 match_base_names=True))
        for (g, a), out in written.items():
            got = json.load(open(out))
            self.assertEqual(self._single(g, a, d + '/ref.json')['index_fp'], got['index_fp'])
            self.assertEqual('precomputed' if '/refseq/' in g else 'hashed',
                             got['files'][0]['digest'])
        # two different digests for one file, one line spelling its path exactly
        g = os.path.join(d, 'refseq', 'graph.dbg')
        conf = os.path.join(d, 'conf.sha256')
        with open(conf, 'w') as out:
            out.write('%s  %s\n%s  refseq/graph.dbg\n' % ('a' * 64, g, 'b' * 64))
        index = im.DigestIndex(im.read_digests([conf]))
        for spelling in (g, os.path.join(d, 'refseq', '.', 'graph.dbg')):
            with self.assertRaises(im.ManifestError):
                index.find(spelling)

    def test_a_pairs_manifest_path_is_kept_on_every_line(self):
        """The manifest path a later line of a pair names is kept, as --write-csv keeps it,
        so that the server accepts the list the tool writes."""
        f = self.files
        csv = os.path.join(self.root, 'in.csv')
        named = os.path.join(self.root, 'b.manifest.json')
        with open(csv, 'w') as out:
            out.write('A,%s,%s\nB,%s,%s,%s\n' % (f['g3.dbg'], f['a3.column.annodbg'],
                                                f['g3.dbg'], f['a3.column.annodbg'], named))
        written = im.batch(_args(server_csv=csv, write_csv=os.path.join(self.root, 'out.csv')))
        self.assertEqual({(f['g3.dbg'], f['a3.column.annodbg']): named}, written)
        self.assertTrue(os.path.isfile(named))
        lines = im.parse_server_csv(os.path.join(self.root, 'out.csv'))
        self.assertEqual([named, named], [c[3] for _, c in lines])
        # two different manifest paths (or index_ns) for one pair are refused
        with open(csv, 'a') as out:
            out.write('C,%s,%s,%s\n' % (f['g3.dbg'], f['a3.column.annodbg'],
                                        os.path.join(self.root, 'c.json')))
        with self.assertRaises(im.ManifestError):
            im.batch(_args(server_csv=csv, force=True))
        with open(csv, 'w') as out:
            out.write('A,%s,%s,,x\nB,%s,%s,,y\n' % (f['g3.dbg'], f['a3.column.annodbg'],
                                                   f['g3.dbg'], f['a3.column.annodbg']))
        with self.assertRaises(im.ManifestError):
            im.batch(_args(server_csv=csv, force=True))

    def test_sidecars_are_those_the_server_loads(self):
        """The bundle is the loader dependency inventory (load_inventory, the mirror of
        metagraph's index_load_inventory): the sidecars the loaders open for the pair, derived
        from the LISTED spelling, and nothing the server does not open. The integration test
        compares it with the binary's `traverse --index-inventory` on real bundles; this one
        pins the rules on stand-ins."""
        d = self.root
        names = ('G.dbg', 'G.dbg.anchors', 'G.dbg.rd_succ', 'G.edgemask', 'G.bloom',
                 'G.dbg.weights', 'X.column.annodbg', 'X.column.annodbg.coords', 'X.seqs',
                 'Y.row_diff_brwt_coord.annodbg', 'Y.seqs', 'Z.row_diff.annodbg',
                 'H.dbg', 'H.bloom')
        for name in names:
            with open(os.path.join(d, name), 'w') as out:
                out.write(name)
        p = lambda n: os.path.join(d, n)
        # a column annotation: no identity sidecar — never the graph's dummy-edge mask or its
        # Bloom filter (derived data: loaded, listed apart by derived_files, in no manifest),
        # the weights, the column's .coords or a .seqs beside an annotation that is not a
        # coordinate one
        self.assertEqual([], im.sidecars(p('G.dbg'), p('X.column.annodbg')))
        self.assertEqual([(p('G.dbg'), 'graph', True),
                          (p('X.column.annodbg'), 'annotation', True)],
                         im.load_inventory(p('G.dbg'), p('X.column.annodbg')))
        # the derived data: the mask when it opens, and only then the Bloom filter
        self.assertEqual([(p('G.edgemask'), 'graph_mask', True, True),
                          (p('G.bloom'), 'graph_bloom', True, True)], im.derived_files(p('G.dbg')))
        self.assertEqual([(p('H.edgemask'), 'graph_mask', False, False),
                          (p('H.bloom'), 'graph_bloom', True, False)], im.derived_files(p('H.dbg')))
        self.assertEqual([], im.derived_files(p('X.column.annodbg')))
        self.assertEqual('graph_mask', im.derived_role('dir/G.edgemask'))
        self.assertEqual('graph_bloom', im.derived_role('other.bloom'))
        self.assertIsNone(im.derived_role('G.dbg'))
        self.assertIsNone(im.derived_role('G.bloom.txt'))
        # a column row-diff annotation: the anchors and fork successors beside the GRAPH,
        # required
        self.assertEqual([(p('G.dbg.anchors'), 'row_diff_anchors', True),
                          (p('G.dbg.rd_succ'), 'row_diff_fork_succ', True)],
                         im.load_inventory(p('G.dbg'), p('Z.row_diff.annodbg'))[-2:])
        # a coordinate annotation: its sequence headers, unless the server runs with
        # --no-coord-mapping (it does not load them then)
        self.assertEqual([p('Y.seqs')],
                         im.sidecars(p('G.dbg'), p('Y.row_diff_brwt_coord.annodbg')))
        self.assertEqual([], im.sidecars(p('G.dbg'), p('Y.row_diff_brwt_coord.annodbg'),
                                         coord_mapping=False))
        self.assertEqual(p('Y.seqs'), im.coordinate_headers(p('Y.row_diff_brwt_coord.annodbg')))
        self.assertIsNone(im.coordinate_headers(p('X.column.annodbg')))
        self.assertEqual([], im.sidecars(p('H.dbg'), p('X.column.annodbg')))
        # required files are listed although missing (the server cannot load without them),
        # and the bundle is then refused, naming them
        self.assertEqual([p('H.dbg.anchors'), p('H.dbg.rd_succ')],
                         im.sidecars(p('H.dbg'), p('Z.row_diff.annodbg')))
        refuse = lambda m: (_ for _ in ()).throw(im.ManifestError(m))
        with self.assertRaises(im.ManifestError) as cm:
            im.bundle_files(p('H.dbg'), p('Z.row_diff.annodbg'), [], fail=refuse)
        self.assertIn('H.dbg.anchors', str(cm.exception))
        # a symlinked spelling: the sidecars beside the links, not beside their targets
        # (two symlinks to one graph and annotation could hide different .seqs files)
        link = os.path.join(d, 'link')
        os.makedirs(link)
        os.symlink(p('G.dbg'), os.path.join(link, 'G.dbg'))
        os.symlink(p('Y.row_diff_brwt_coord.annodbg'),
                   os.path.join(link, 'Y.row_diff_brwt_coord.annodbg'))
        self.assertEqual([], im.sidecars(os.path.join(link, 'G.dbg'),
                                         os.path.join(link, 'Y.row_diff_brwt_coord.annodbg')))
        with open(os.path.join(link, 'Y.seqs'), 'w') as out:
            out.write('other headers')
        self.assertEqual([os.path.join(link, 'Y.seqs')],
                         im.sidecars(os.path.join(link, 'G.dbg'),
                                     os.path.join(link, 'Y.row_diff_brwt_coord.annodbg')))
        # a further graph or annotation is not part of one index's bundle
        with self.assertRaises(im.ManifestError):
            im.bundle_files(p('G.dbg'), p('X.column.annodbg'),
                            [p('Y.row_diff_brwt_coord.annodbg')], fail=refuse)
        # nor the graph's derived data, as an extra file
        for extra in ('G.edgemask', 'G.bloom', 'H.bloom'):
            with self.assertRaises(im.ManifestError) as cm:
                im.bundle_files(p('G.dbg'), p('X.column.annodbg'), [p(extra)], fail=refuse)
            self.assertIn(im.DERIVED_RULE, str(cm.exception))

    def test_derived_data_is_never_in_a_manifest(self):
        """The graph's mask and Bloom filter are derived data, not part of index_fp. A
        manifest written with them beside the graph is the one written without them (the
        same index_fp: adding a mask to a deployed index leaves its identity unchanged), in
        single and in batch mode; --verify flags a manifest that lists one (the server
        refuses it)."""
        f = self.files
        g, a = f['g1.dbg'], f['a2.column.annodbg']
        with_derived = self._single(g, a, os.path.join(self.root, 'with.json'))
        self.assertEqual(['a2.column.annodbg', 'g1.dbg'],
                         [e['path'] for e in with_derived['files']])
        for name in ('g1.edgemask', 'g1.bloom'):
            os.rename(f[name], f[name] + '.away')
        without = self._single(g, a, os.path.join(self.root, 'without.json'))
        for name in ('g1.edgemask', 'g1.bloom'):
            os.rename(f[name] + '.away', f[name])
        self.assertEqual(without, with_derived)
        csv = os.path.join(self.root, 'one.csv')
        with open(csv, 'w') as out:
            out.write('B,%s,%s\n' % (g, a))
        written = im.batch(_args(server_csv=csv, out_dir=os.path.join(self.root, 'b')))
        with open(written[(g, a)]) as fh:
            self.assertEqual(without['index_fp'], json.load(fh)['index_fp'])
        # --verify: the manifest as written passes; with the mask listed it is flagged
        manifest = os.path.join(self.root, 'with.json')
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            im.verify(manifest, root=self.root)
        self.assertIn('OK', out.getvalue())
        listed = dict(with_derived)
        with open(f['g1.edgemask'], 'rb') as fh:
            data = fh.read()
        entries = list(with_derived['files']) + [{'path': 'g1.edgemask', 'size': len(data),
                                                  'sha256': hashlib.sha256(data).hexdigest()}]
        entries.sort(key=lambda e: e['path'].encode())
        listed['files'] = entries
        listed['index_fp'] = hashlib.sha256(im.canonical_lines(entries).encode()).hexdigest()
        bad = os.path.join(self.root, 'listed.json')
        with open(bad, 'w') as fh:
            json.dump(listed, fh)
        out = io.StringIO()
        with contextlib.redirect_stdout(out), self.assertRaises(SystemExit):
            im.verify(bad, root=self.root)
        self.assertIn('derived data listed: g1.edgemask (' + im.DERIVED_RULE, out.getvalue())
        # (the digests and the stated index_fp are right: only the rule is broken)
        self.assertEqual(1, out.getvalue().count('\n'))

    def test_no_coord_mapping_leaves_the_headers_out_in_both_modes(self):
        """--no-coord-mapping: the server does not load the .seqs, so neither single nor batch
        mode lists it, and the two write the same manifest."""
        f = self.files
        g, a = f['g1.dbg'], f['a1.relaxed.row_diff_brwt_coord.annodbg']
        single = os.path.join(self.root, 'single.json')
        im.write(_args(graph=g, annotation=a, output=single, no_coord_mapping=True))
        with open(single) as fh:
            ref = json.load(fh)
        self.assertNotIn('a1.relaxed.seqs', [e['path'] for e in ref['files']])
        csv = os.path.join(self.root, 'one.csv')
        with open(csv, 'w') as out:
            out.write('A,%s,%s\n' % (g, a))
        written = im.batch(_args(server_csv=csv, out_dir=self.root + '/nc',
                                 no_coord_mapping=True))
        with open(written[(g, a)]) as fh:
            got = json.load(fh)
        self.assertEqual(ref['index_fp'], got['index_fp'])
        with_headers = os.path.join(self.root, 'with.json')
        im.write(_args(graph=g, annotation=a, output=with_headers))
        with open(with_headers) as fh:
            self.assertIn('a1.relaxed.seqs', [e['path'] for e in json.load(fh)['files']])

    def test_lines_the_server_refuses_are_refused(self):
        for line in ('A,g,a,m,ns,extra', 'A,g', 'A,g,a,m,bad ns'):
            path = os.path.join(self.root, 'bad.csv')
            with open(path, 'w') as out:
                out.write(line + '\n')
            with self.assertRaises(im.ManifestError):
                im.parse_server_csv(path)
        # two pairs writing one manifest
        path = os.path.join(self.root, 'clash.csv')
        f = self.files
        with open(path, 'w') as out:
            out.write('A,%s,%s,%s/m.json\nB,%s,%s,%s/m.json\n'
                      % (f['g1.dbg'], f['a2.column.annodbg'], self.root,
                         f['g3.dbg'], f['a3.column.annodbg'], self.root))
        with self.assertRaises(im.ManifestError):
            im.batch(_args(server_csv=path))
        # a missing bundle file names its line
        path = os.path.join(self.root, 'missing.csv')
        with open(path, 'w') as out:
            out.write('A,%s/nope.dbg,%s\n' % (self.root, f['a2.column.annodbg']))
        with self.assertRaises(im.ManifestError) as cm:
            im.batch(_args(server_csv=path, out_dir=self.root + '/o'))
        self.assertIn('line 1', str(cm.exception))


if __name__ == '__main__':
    unittest.main()
