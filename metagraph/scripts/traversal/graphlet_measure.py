#!/usr/bin/env python3
"""Measure MGT graphlets (DESIGN-traverse-graphlet.md §7): the size estimates of the
design are estimates under an interning assumption, not bounds, and are to be replaced
by these numbers on real retrievals (the SRA case).

Per document: bytes and lines per record type, gzip size, parse / dump / to_json time,
the model's memory estimate, label-set cardinalities (G entry/end, P), how many sets
are distinct (the interner's hit rate), and membership changes per step (|delta|
between consecutive P sets and between a parent's end set and a child's entry set).
With --full the size of the detail: full JSON is reported beside the body.

Inputs: .mgt files, detail: graphlet response JSON files (every result's graphlet),
--server HOST:PORT --request REQUEST.json, or --cli METAGRAPH --graph G --annotation A
--request REQUEST.json (each run once with detail graphlet and once with detail full;
--request may be given several times).

usage:
    graphlet_measure.py FILE [FILE ...] [--json]
    graphlet_measure.py --server localhost:5555 --request req.json [--json]
    graphlet_measure.py --cli build/metagraph --graph G.dbg --annotation A.annodbg \
        [--index-manifest M.json] --request req.json [--request ...] [--json]
"""

import argparse
import copy
import gzip
import json
import os
import statistics
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(HERE)), 'api', 'python'))

from metagraph.traverse import from_response  # noqa: E402
from metagraph.traverse import parser as _parser  # noqa: E402


def _stats(xs):
    if not xs:
        return {'n': 0}
    return {'n': len(xs), 'mean': round(statistics.fmean(xs), 2), 'max': max(xs),
            'median': statistics.median(xs)}


def measure_text(text, full_json=None, envelope=None):
    data = text.encode('utf-8')
    per = {}
    for line in text.split('\n')[:-1]:
        tag = line[:1]
        b = per.setdefault(tag, [0, 0])
        b[0] += 1
        b[1] += len(line.encode('utf-8')) + 1
    t0 = time.perf_counter()
    reader = _parser._Reader(text.split('\n')[:-1])
    g = reader.run()
    t1 = time.perf_counter()
    dumped = g.dump()
    t2 = time.perf_counter()
    out = {
        'bytes': len(data), 'gzip_bytes': len(gzip.compress(data, 9)),
        'lines': text.count('\n'),
        'records': {k: {'count': v[0], 'bytes': v[1]} for k, v in sorted(per.items())},
        'parse_ms': round((t1 - t0) * 1000, 2), 'dump_ms': round((t2 - t1) * 1000, 2),
        'canonical': dumped == text, 'memory_bytes': g.memory_bytes(),
        'interner': {'distinct_sets': len(reader.interner), 'hits': reader.interner.hits,
                     'misses': reader.interner.misses,
                     'hit_rate': round(reader.interner.hits
                                       / max(1, reader.interner.hits + reader.interner.misses), 3),
                     'ids_stored': reader.interner.ids_stored},
        'arms': {},
    }
    for side, a in g.arms.items():
        entry = [len(s.entry) for s in a.segments]
        end = [len(s.end) for s in a.segments]
        pres = [len(p.labels) for s in a.segments for p in s.presence]
        p_delta, edge_delta = [], []
        for s in a.segments:
            prev = set(s.entry)
            for p in s.presence:
                cur = set(p.labels)
                p_delta.append(len(prev ^ cur))
                prev = cur
            for c in s.children:
                edge_delta.append(len(set(s.end) ^ set(a.segments[c].entry)))
        steps = sum(s.length_bp for s in a.segments) or 1
        out['arms'][side] = {
            'segments': len(a.segments), 'runs': len(a.runs),
            'leaves': a.counts.leaves, 'splits': a.counts.splits, 'merges': a.counts.merges,
            'bases': a.counts.bases, 'entry_set': _stats(entry), 'end_set': _stats(end),
            'presence_set': _stats(pres),
            'membership_changes_per_step': round((sum(p_delta) + sum(edge_delta)) / steps, 4),
            'p_delta': _stats(p_delta), 'edge_delta': _stats(edge_delta)}
    if envelope is not None:
        # |envelope| is the graphlet built with its envelope (from_response)
        t3 = time.perf_counter()
        j = envelope.to_json()
        out['to_json_ms'] = round((time.perf_counter() - t3) * 1000, 2)
        out['to_json_compact_bytes'] = len(json.dumps(j, separators=(',', ':')).encode())
    if full_json is not None:
        compact = json.dumps(full_json, separators=(',', ':')).encode('utf-8')
        out['full_json_bytes'] = len(compact)
        out['full_json_gzip_bytes'] = len(gzip.compress(compact, 9))
        out['body_vs_full'] = round(len(data) / max(1, len(compact)), 4)
    return out


def measure_response(resp, full=None):
    outs = []
    for i, r in enumerate(resp.get('results', [])):
        if 'graphlet' not in r:
            outs.append({'seed': i, 'error': r.get('error')})
            continue
        g = from_response(r, resp)
        full_result = full['results'][i] if full else None
        m = measure_text(r['graphlet'], full_result, g)
        m['seed'] = i
        m['summary_bytes'] = len(json.dumps({k: v for k, v in r.items() if k != 'graphlet'},
                                            separators=(',', ':')).encode())
        outs.append(m)
    return outs


def run_cli(args, request):
    """`metagraph traverse` on one request -> the parsed response."""
    with tempfile.NamedTemporaryFile('w', suffix='.json', delete=False) as f:
        json.dump(request, f)
    try:
        cmd = [args.cli, 'traverse', '-i', args.graph, '-a', args.annotation]
        if args.index_manifest:
            cmd += ['--index-manifest', args.index_manifest]
        p = subprocess.run(cmd + [f.name], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    finally:
        os.unlink(f.name)
    resp = json.loads(p.stdout.decode('utf-8'))
    if p.returncode:
        raise SystemExit('traverse failed: %s' % resp.get('error'))
    return resp


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('files', nargs='*')
    ap.add_argument('--server')
    ap.add_argument('--cli', metavar='METAGRAPH')
    ap.add_argument('--graph')
    ap.add_argument('--annotation')
    ap.add_argument('--index-manifest')
    ap.add_argument('--request', action='append', default=[])
    ap.add_argument('--json', action='store_true', help='machine-readable output')
    args = ap.parse_args(argv)
    results = []
    for path in args.request if (args.server or args.cli) else []:
        with open(path) as f:
            req = json.load(f)
        resp = {}
        for detail in ('graphlet', 'full'):
            r = copy.deepcopy(req)
            r.setdefault('strategy', {}).setdefault('output', {})['detail'] = detail
            if args.server:
                from metagraph.traverse import TraverseClient
                host, port = args.server.rsplit(':', 1)
                resp[detail] = TraverseClient(host, int(port)).traverse_raw(r)
            else:
                resp[detail] = run_cli(args, r)
        results.append({'source': path, 'seeds': measure_response(resp['graphlet'],
                                                                  resp['full'])})
    for path in args.files:
        with open(path, 'r', encoding='utf-8', newline='') as f:
            text = f.read()
        if path.endswith('.json'):
            results.append({'source': path, 'seeds': measure_response(json.loads(text))})
        else:
            results.append({'source': path, 'seeds': [measure_text(text)]})
    if args.json:
        json.dump(results, sys.stdout, indent=1)
        print()
        return 0
    for r in results:
        for m in r['seeds']:
            if 'error' in m:
                print('%s: seed %s failed: %s' % (r['source'], m['seed'], m['error']))
                continue
            print('%s: %d bytes (%d gzip), %d lines, parse %.1f ms, model ~%d KB, '
                  'interner %d sets (hit rate %.2f)'
                  % (r['source'], m['bytes'], m['gzip_bytes'], m['lines'], m['parse_ms'],
                     m['memory_bytes'] // 1024, m['interner']['distinct_sets'],
                     m['interner']['hit_rate']))
            for tag, v in m['records'].items():
                print('   %s %6d records %9d bytes' % (tag, v['count'], v['bytes']))
            for side, a in m['arms'].items():
                print('   %s: %d segments, %d runs, %d leaves; entry sets %s; '
                      'membership changes/step %s'
                      % (side, a['segments'], a['runs'], a['leaves'], a['entry_set'],
                         a['membership_changes_per_step']))
            if 'full_json_bytes' in m:
                print('   detail full: %d bytes compact (%d gzip): body/full = %.3f'
                      % (m['full_json_bytes'], m['full_json_gzip_bytes'], m['body_vs_full']))
    return 0


if __name__ == '__main__':
    sys.exit(main())
