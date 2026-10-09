"""The record-coordinate snapshots of the library's tests: real server responses on
mini_refseq with strategy.output.coordinates (feature level 6, DESIGN §18), committed
under api/python/tests/data/traverse/coords/ so that test_traverse_coordinates.py and the
golden file golden/coordinates.json.gz run with no server.

Each cell is <name>.request.json.gz (the request, detail left out), <name>.graphlet.json.gz
and <name>.full.json.gz (the same request in detail graphlet and full), gzipped with mtime 0
so that a rebuild with the same binary is byte-identical. The cells cover what the block
can state: cut lists (cap 1) with the seed-level limitation and its K record, "unlimited",
a run marked lower_bound (a switch into a live label), the column and mixed kinds, both
arms and one, the null reasons "support kmer" and "partial derivation" (D3), a D3 result
without coordinates, and a memory stop offering drop_coordinates.

    python3 scripts/traversal/traverse_coords_snapshots.py --binary BUILD/metagraph [--out DIR]
                                                           [--index DIR]

The index is build/mini_refseq (scripts/traversal/build_mini_refseq.sh); the server runs
with --traverse-chunk-target-ms 0, so the D3 cells (a time budget of 1e-6 ms) stop after the
first k-mer deterministically. Rebuilding changes the snapshots only where the server's
output changed: then record the golden file again (api/python/tests/traverse_golden.py --out
data/traverse/golden/coordinates.json.gz data/traverse/coords) and say why.
"""

import argparse
import copy
import gzip
import json
import os
import socket
import subprocess
import sys
import time
import urllib.error
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
OUT = os.path.join(REPO, 'api', 'python', 'tests', 'data', 'traverse', 'coords')

# seeds of mini_refseq: a 150 bp window of NZ_CP030345.1 that the index holds six times
# (a repeat, so lists can be cut), blaNDM-1 (813 bp), and a 400 bp seed whose derivation
# the D3 cells interrupt
REP0 = ('CAAAGTTAGCGATGAGGCAGCCTTTTGTCTTATTCAAAGGCCTTACATTTCAAAAACTCTGCTTACCAGGCGCATTTCGC'
        'CCAGGGGATCACCATAATAAAATGCTGAGGCCTGGCCTTTGCGTAGTGCACGCATCACCTCAATACCTTT')
D3 = ('GTAAAGCCCTCGCTAGATTTTAATGCGGATGTTGCGATTACTTCGCCAACTATTGCGATAACAAGAAAAAGCCAGCCTTT'
      'CATGATATATCTCCCAATTTGTGTAGGGCTTATTATGCACGCTTGCTTACACCATTAGAGAAATTTGCTCAGCTTGTTGAT'
      'TATCATATGGCTTTTGAAACTGTCGCACCTCATGTTTGAATTCGCCCCATATTTTTGCTACAGTGAACCAAATTAAGATCAT'
      'CTATTTACTAGGCCTCGCATTTGCGGGGTTTTTAATGCTGAATAAAAGGAAAACTTGATGGAATTGCCCAATATTATGCAC'
      'CCGGTCGCGAAGCTGAGCACCGCATTAGCCGCTGCATTGATGCTGAGCGGGTGCATGCCCGGTGAAATCCGCCCGA')

_TRACE = {'support': 'trace', 'branching': {'on_reconverge': 'keep', 'max_label_branches': 2}}
_SWITCH = {'support': 'trace',
           'labels': {'change_cost': {'model': 'constant', 'value': 0.5}, 'loss_budget': 2},
           'branching': {'on_reconverge': 'keep', 'max_label_branches': 2}}


def _strategy(base, ext=300, direction=None, bounds=None, labels=None, **output):
    s = copy.deepcopy(base)
    s['bounds'] = dict({'max_extension_bp': ext}, **(bounds or {}))
    if direction:
        s['direction'] = direction
    if labels:
        s.setdefault('labels', {}).update(labels)
    s['output'] = dict({'timing': False}, **output)
    return s


# name -> (seeds, strategy, what it covers)
CELLS = {
    'rep0_cut1': ([{'seed_id': 'rep0', 'sequence': REP0}],
                  _strategy(_TRACE, coordinates=True, max_coordinate_occurrences=1),
                  'both arms, header labels, cap 1: 15 lists cut, the K record'),
    'rep0_unlimited_left': ([{'seed_id': 'rep0', 'sequence': REP0}],
                            _strategy(_TRACE, direction='left', coordinates=True,
                                      max_coordinate_occurrences='unlimited'),
                            'the left arm only, "unlimited": nothing cut'),
    'rep0_lower_bound': ([{'seed_id': 'rep0', 'sequence': REP0}],
                         _strategy(_SWITCH, direction='both', coordinates=True),
                         'switches (cost 0.5, budget 2): a run marked lower_bound'),
    'rep0_column': ([{'seed_id': 'rep0', 'sequence': REP0}],
                    _strategy(_SWITCH, 600, 'both', labels={'seed_label_kind': 'column'},
                              coordinates=True),
                    'column labels with switches: kind column, trace_record_boundaries'),
    'rep0_mixed': ([{'seed_id': 'rep0', 'sequence': REP0}],
                   _strategy(_SWITCH, labels={'extra': ['562']}, coordinates=True),
                   'header seed labels and an extra column label: kind mixed'),
    'rep0_kmer_null': ([{'seed_id': 'rep0', 'sequence': REP0}],
                       _strategy({}, 100, coordinates=True),
                       'support kmer: null, reason "support kmer"'),
    'd3_trace_coordinates': ([{'sequence': D3}],
                             _strategy(_TRACE, 100, bounds={'time_budget_ms': 1e-6},
                                       coordinates=True),
                             'D3 under trace: null, reason "partial derivation"'),
    'd3_kmer': ([{'sequence': D3}],
                _strategy({}, 100, bounds={'time_budget_ms': 1e-6}),
                'D3 without coordinates: a walked result with a derivation limitation'),
    'd3_memory_stop': ([{'sequence': D3}],
                       _strategy(_TRACE, 1500, bounds={'max_memory_mb': 1},
                                 coordinates=True),
                       'a memory stop offering drop_coordinates'),
}


def _free_port():
    s = socket.socket()
    s.bind(('127.0.0.1', 0))
    port = s.getsockname()[1]
    s.close()
    return port


class _Server:
    def __init__(self, binary, index):
        port = _free_port()
        cmd = [binary, 'server_query', '-i', os.path.join(index, 'graph_k31.dbg'),
               '-a', os.path.join(index, 'annotation.relaxed.relabeled.row_diff_brwt_coord'
                                         '.annodbg'),
               '--port', str(port), '--address', '127.0.0.1', '-p', '2', '--mmap',
               '--traverse-max-time-ms', '120000', '--traverse-chunk-target-ms', '0',
               '--index-name', 'mini_refseq', '--index-manifest',
               os.path.join(index, 'annotation.relaxed.relabeled.manifest.json')]
        self.p = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        self.url = 'http://127.0.0.1:%d' % port
        for _ in range(600):
            try:
                with urllib.request.urlopen(self.url + '/traverse/capabilities', timeout=5):
                    return
            except (OSError, urllib.error.URLError):
                if self.p.poll() is not None:
                    raise RuntimeError('the server exited')
                time.sleep(0.2)
        raise RuntimeError('the server did not start')

    def post(self, payload):
        req = urllib.request.Request(self.url + '/traverse', json.dumps(payload).encode(),
                                     {'Content-Type': 'application/json'})
        try:
            with urllib.request.urlopen(req, timeout=600) as r:
                return r.status, r.read()
        except urllib.error.HTTPError as e:
            return e.code, e.read()

    def close(self):
        self.p.kill()
        self.p.wait()


def _write_gz(path, data):
    with open(path, 'wb') as f:
        with gzip.GzipFile(fileobj=f, mode='wb', mtime=0, filename='') as z:
            z.write(data)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--binary', required=True)
    ap.add_argument('--index', default=os.path.join(REPO, 'build', 'mini_refseq'))
    ap.add_argument('--out', default=OUT)
    args = ap.parse_args(argv)
    os.makedirs(args.out, exist_ok=True)
    s = _Server(args.binary, args.index)
    try:
        for name, (seeds, strategy, _) in CELLS.items():
            _write_gz(os.path.join(args.out, name + '.request.json.gz'),
                      json.dumps({'seeds': seeds, 'strategy': strategy}, sort_keys=True,
                                 indent=1).encode())
            for detail in ('graphlet', 'full'):
                st = copy.deepcopy(strategy)
                st['output']['detail'] = detail
                code, body = s.post({'seeds': seeds, 'strategy': st})
                if code != 200:
                    raise SystemExit('%s %s: HTTP %d %s' % (name, detail, code, body[:300]))
                _write_gz(os.path.join(args.out, '%s.%s.json.gz' % (name, detail)), body)
            print(name, 'ok')
    finally:
        s.close()


if __name__ == '__main__':
    sys.exit(main())
