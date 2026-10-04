"""parser.standalone_text (pass 5, W7, R14): the standalone .mgt text of one seed — H, the J
envelope and the body — spliced without parsing the body, byte for byte what
dump(from_response(result, response), envelope=True) writes; the J envelope's usage reduced to
the totals plus this seed's per_seed entry, the same rule in save(), dump() and the store; the
envelope names `view` and `derived_from` reserved."""
import copy
import glob
import gzip
import json
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import traverse_testlib  # noqa: E402,F401  (puts the library on the path)
from metagraph.traverse import parser  # noqa: E402
from metagraph.traverse.parser import (GraphletFormatError, dump, from_response, parse,  # noqa: E402
                                       save, seed_envelope, standalone_text)
from metagraph.traverse.store import GraphletStore  # noqa: E402

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'traverse')


def _responses():
    files = sorted(glob.glob(ROOT + '/documents/**/*.graphlet.json', recursive=True)) + \
        sorted(glob.glob(ROOT + '/real/*/*.graphlet.json.gz'))
    for f in files:
        op = gzip.open if f.endswith('.gz') else open
        with op(f, 'rt', encoding='utf-8') as fh:
            yield f, json.load(fh)


def _with_usage(n=3):
    """A response of |n| seeds with an attempt's usage, built from a fixture's first result
    (each copy's H record names its own seed_index)."""
    with open(ROOT + '/documents/fork.graphlet.json', encoding='utf-8') as f:
        resp = json.load(f)
    r0 = resp['results'][0]
    resp['results'] = []
    for i in range(n):
        r = copy.deepcopy(r0)
        head, rest = r['graphlet'].split('\n', 1)
        t = head.split(' ')
        t[11] = str(i)
        r['graphlet'] = ' '.join(t) + '\n' + rest
        resp['results'].append(r)
    resp['usage'] = {'attempt_id': 'a-17', 'budget_id': 'b-3', 'server_instance': 'f' * 16,
                     'reason': 'completed', 'work_units': 3 * 7,
                     'seeds': {'requested': n, 'started': n, 'finished': n, 'abandoned': 0},
                     'per_seed': [{'index': i, 'outcome': 'complete', 'work_units': 7}
                                  for i in range(n)]}
    return resp


class TestStandaloneText(unittest.TestCase):
    def test_byte_equal_to_dump_on_every_fixture(self):
        n = 0
        for f, resp in _responses():
            for i, r in enumerate(resp.get('results', [])):
                if 'graphlet' not in r:
                    continue
                n += 1
                got = standalone_text(r['graphlet'], r, resp)
                ref = dump(from_response(r, resp), envelope=True)
                self.assertEqual(ref, got, '%s result %d' % (f, i))
                # bytes, as a body may arrive
                self.assertEqual(got, standalone_text(r['graphlet'].encode('utf-8'), r, resp))
        self.assertGreaterEqual(n, 34)

    def test_usage_is_reduced_to_this_seed(self):
        resp = _with_usage()
        for i, r in enumerate(resp['results']):
            text = standalone_text(r['graphlet'], r, resp)
            self.assertEqual(dump(from_response(r, resp), envelope=True), text)
            j = json.loads(text.split('\n')[1][2:])
            usage = j['usage']
            self.assertEqual([{'index': i, 'outcome': 'complete', 'work_units': 7}],
                             usage['per_seed'])
            # the totals are kept
            self.assertEqual(21, usage['work_units'])
            self.assertEqual('a-17', usage['attempt_id'])
            # the response itself is not changed
            self.assertEqual(3, len(resp['usage']['per_seed']))
            # canonical, and idempotent through parse and dump
            self.assertEqual(text, dump(parse(text)))
            g = parse(text)
            self.assertEqual([{'index': i, 'outcome': 'complete', 'work_units': 7}],
                             g.envelope['usage']['per_seed'])
        # seed_envelope is idempotent and leaves an envelope without per_seed as it is
        env = {k: v for k, v in resp.items() if k != 'results'}
        once = seed_envelope(env, 1)
        self.assertEqual(once, seed_envelope(once, 1))
        self.assertEqual({'a': 1}, seed_envelope({'a': 1}, 0))
        self.assertEqual({}, seed_envelope(None, 0))

    def test_save_and_the_store_keep_the_same_bytes(self):
        resp = _with_usage()
        with tempfile.TemporaryDirectory() as d:
            store = GraphletStore(os.path.join(d, 'store'))
            for i, r in enumerate(resp['results']):
                text = standalone_text(r['graphlet'], r, resp)
                path = os.path.join(d, 'seed%d.mgt' % i)
                save(from_response(r, resp), path)
                with open(path, encoding='utf-8', newline='') as f:
                    self.assertEqual(text, f.read())
                h = store.put(resp, None, i)
                # the stored entry keeps this seed's usage only
                entry = store.get(h)
                self.assertEqual([i], [e['index'] for e in entry.envelope['usage']['per_seed']])
                out = os.path.join(d, 'stored%d.mgt' % i)
                store.save(h, out)
                with open(out, encoding='utf-8', newline='') as f:
                    self.assertEqual(text, f.read())
                # read back from disk (not the resident model): the same text
                store._drop_ram(h)
                store.save(h, out)
                with open(out, encoding='utf-8', newline='') as f:
                    self.assertEqual(text, f.read())

    def test_refusals(self):
        resp = _with_usage(1)
        r = resp['results'][0]
        for name in parser.RESERVED_ENVELOPE_NAMES:
            bad = dict(resp, **{name: {'x': 1}})
            with self.assertRaises(ValueError):
                standalone_text(r['graphlet'], r, bad)
            with self.assertRaises(ValueError):
                from_response(r, bad)
        # cut in transport: bytes and lines, as from_response
        cut = dict(r, graphlet=r['graphlet'][:-10])
        with self.assertRaises(GraphletFormatError):
            standalone_text(cut['graphlet'], cut, resp)
        with self.assertRaises(GraphletFormatError):
            from_response(cut, resp)
        lines = dict(r, graphlet_lines=r['graphlet_lines'] + 1)
        with self.assertRaises(GraphletFormatError):
            standalone_text(r['graphlet'], lines, resp)
        # a body that already carries a J line is a saved file, not a server body
        saved = standalone_text(r['graphlet'], r, resp)
        nocount = {k: v for k, v in r.items() if k not in ('graphlet_bytes', 'graphlet_lines')}
        with self.assertRaises(ValueError):
            standalone_text(saved, nocount, resp)
        # no H, no Z, no final line feed
        body = r['graphlet']
        for broken in ('X' + body[1:], body[:body.rfind('Z ')] + 'Q 3\n', body[:-1]):
            with self.assertRaises(GraphletFormatError):
                standalone_text(broken, nocount, resp)

    def test_a_file_saved_before_the_reduced_envelope_still_loads(self):
        """Review of pass 5: a file the stage-4 library saved kept every seed's per_seed entry
        in its J line. It still parses and loads; it is not canonical under the reduced rule
        (dump(parse(text)) differs from it), and save() rewrites it reduced (SPEC §7.5.2)."""
        from metagraph.traverse.parser import is_canonical, load
        resp = _with_usage()
        r = resp['results'][1]
        new = standalone_text(r['graphlet'], r, resp)
        lines = new.split('\n')
        j = json.loads(lines[1][2:])
        j['usage'] = copy.deepcopy(resp['usage'])          # as the stage-4 library wrote it
        old = '\n'.join([lines[0], 'J ' + json.dumps(j, sort_keys=True, separators=(',', ':'),
                                                       ensure_ascii=False)] + lines[2:])
        self.assertNotEqual(new, old)
        g = parse(old)
        self.assertEqual(3, len(g.envelope['usage']['per_seed']))
        self.assertFalse(is_canonical(old))
        self.assertTrue(is_canonical(new))
        self.assertEqual(new, dump(g, envelope=True))
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, 'old.mgt')
            with open(path, 'w', encoding='utf-8') as f:
                f.write(old)
            loaded = load(path)
            save(loaded, os.path.join(d, 'again.mgt'))
            with open(os.path.join(d, 'again.mgt'), encoding='utf-8') as f:
                self.assertEqual(new, f.read())


if __name__ == '__main__':
    unittest.main()
