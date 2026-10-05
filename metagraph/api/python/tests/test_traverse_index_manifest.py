"""scripts/traversal/index_manifest.py, batch mode (pass 5, W2): one manifest per (graph,
annotation) pair of a multi-graph server's list, each distinct file hashed once by parallel
streams or taken from precomputed sha256 digests, equal to what single mode writes."""
import argparse
import contextlib
import hashlib
import io
import importlib.util
import json
import os
import tempfile
import threading
import unittest

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
             name=None, extra=[], no_coord_mapping=False, inventory=False)
    a.update(kw)
    return argparse.Namespace(**a)


class TestIndexManifestBatch(unittest.TestCase):
    def setUp(self):
        # the script reports to stdout and stderr; the tests read what it writes
        self._quiet = contextlib.ExitStack()
        self._quiet.enter_context(contextlib.redirect_stdout(io.StringIO()))
        self._quiet.enter_context(contextlib.redirect_stderr(io.StringIO()))
        self.addCleanup(self._quiet.close)
        self.tmp = tempfile.TemporaryDirectory()
        d = self.root = self.tmp.name
        # two bundles sharing one graph (one graph, two annotations), and a third graph. The
        # shared graph's mask and Bloom filter are loaded for both of its pairs, the coordinate
        # annotation's .seqs for its own; the files the server never opens for these pairs
        # (row-diff anchors without a .row_diff annotation, weights, a column's .coords) are
        # in no bundle and never hashed
        self.files = {}
        self.unloaded = {'g1.dbg.anchors', 'g1.dbg.weights', 'a2.column.annodbg.coords'}
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

    def _single(self, graph, annotation, out):
        im.write(_args(graph=graph, annotation=annotation, output=out))
        with open(out) as f:
            return json.load(f)

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
        """Review of pass 5: bundles share base names (graph.dbg). A checksum file written in
        one bundle's directory with relative paths named the other bundle's files by base name,
        and their manifest took the wrong digests; a line spelling a path exactly also skipped
        the check for a second, different digest of the same file."""
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
        """Review of pass 5: the manifest path a later line of a pair named was ignored, while
        --write-csv kept it, and the server refused the list the tool wrote."""
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
        metagraph's index_load_inventory; review of pass 5, findings 2 and 3): the sidecars
        the loaders open for the pair, derived from the LISTED spelling, and nothing the
        server does not open. The integration test compares it with the binary's
        `traverse --index-inventory` on real bundles; this one pins the rules on stand-ins."""
        d = self.root
        names = ('G.dbg', 'G.dbg.anchors', 'G.dbg.rd_succ', 'G.edgemask', 'G.bloom',
                 'G.dbg.weights', 'X.column.annodbg', 'X.column.annodbg.coords', 'X.seqs',
                 'Y.row_diff_brwt_coord.annodbg', 'Y.seqs', 'Z.row_diff.annodbg',
                 'H.dbg', 'H.bloom')
        for name in names:
            with open(os.path.join(d, name), 'w') as out:
                out.write(name)
        p = lambda n: os.path.join(d, n)
        # a column annotation: the graph's dummy-edge mask and, as the mask opens, its Bloom
        # filter (finding 3: loaded, and able to change answers); never the weights, the
        # column's .coords or a .seqs beside an annotation that is not a coordinate one
        self.assertEqual([p('G.edgemask'), p('G.bloom')],
                         im.sidecars(p('G.dbg'), p('X.column.annodbg')))
        self.assertEqual([(p('G.dbg'), 'graph', True), (p('G.edgemask'), 'graph_mask', False),
                          (p('G.bloom'), 'graph_bloom', False),
                          (p('X.column.annodbg'), 'annotation', True)],
                         im.load_inventory(p('G.dbg'), p('X.column.annodbg')))
        # a column row-diff annotation: the anchors and fork successors beside the GRAPH,
        # required
        self.assertEqual([(p('G.dbg.anchors'), 'row_diff_anchors', True),
                          (p('G.dbg.rd_succ'), 'row_diff_fork_succ', True)],
                         im.load_inventory(p('G.dbg'), p('Z.row_diff.annodbg'))[-2:])
        # a coordinate annotation: its sequence headers, unless the server runs with
        # --no-coord-mapping (it does not load them then)
        self.assertEqual([p('G.edgemask'), p('G.bloom'), p('Y.seqs')],
                         im.sidecars(p('G.dbg'), p('Y.row_diff_brwt_coord.annodbg')))
        self.assertEqual([p('G.edgemask'), p('G.bloom')],
                         im.sidecars(p('G.dbg'), p('Y.row_diff_brwt_coord.annodbg'),
                                     coord_mapping=False))
        self.assertEqual(p('Y.seqs'), im.coordinate_headers(p('Y.row_diff_brwt_coord.annodbg')))
        self.assertIsNone(im.coordinate_headers(p('X.column.annodbg')))
        # the Bloom filter only after the mask: DBGSuccinct::load reads it only then
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
        # (finding 2: two symlinks to one graph and annotation hid different .seqs files)
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
