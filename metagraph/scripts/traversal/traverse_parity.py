#!/usr/bin/env python3
"""
Parity of the /traverse and /resolve routes with the labeled-traversal prototype
(gr/labeled-traversal, d5d26baa) on the v1 feature set.

The requests (REQUESTS below) are v1 requests on the mini RefSeq index
(fetch_mini_refseq.sh). Both binaries' responses are projected onto the fields v1 states;
the prototype's projections are committed as golden files (tests/data/traverse_parity/).

    # check this build against the golden files (CI): starts a server on the index
    traverse_parity.py --metagraph build/metagraph --index build/mini_refseq

    # record the golden files with a build of the prototype, through its `traverse` command
    traverse_parity.py --metagraph <d5d26baa>/build/metagraph --index build/mini_refseq --record

The prototype follows a header label as it follows a column label: by membership. v1 follows a
header label along its record's coordinates. So the requests are column walks, and header
walks only where no record ends or repeats within the radius. Header walks are checked
against the source records in tests/graph/traversal/test_mini_refseq.cpp.

The prototype is sent each request translated into its schema with the knobs v1 fixes stated
explicitly: no reconvergence merging (paths are kept apart), hairpins skipped, 64 splits per
path, no continuation sequences, detail full.
"""
import argparse
import json
import os
import socket
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
GOLDEN_DIR = os.path.join(HERE, '..', '..', 'tests', 'data', 'traverse_parity')
GRAPH = 'graph_k31.dbg'
ANNOTATION = 'annotation.relaxed.relabeled.row_diff_brwt_coord.annodbg'
TIME_BUDGET_MS = 600_000


def read_records(index):
    records = {}
    fasta = os.path.join(index, 'fasta')
    for name in sorted(os.listdir(fasta)):
        header = None
        with open(os.path.join(fasta, name)) as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    header = line[1:].split()[0]
                    records[header] = []
                elif header:
                    records[header].append(line)
    return {header: ''.join(lines) for header, lines in records.items()}


def requests(index):
    """name -> (kind, v1 request)"""
    with open(os.path.join(index, 'sanity', 'ndm1.fa')) as f:
        ndm1 = ''.join(line.strip() for line in f if not line.startswith('>'))
    records = read_records(index)
    names = sorted(records)

    def traverse(seed, **strategy):
        strategy.setdefault('bounds', {})
        strategy['bounds'].setdefault('time_budget_ms', TIME_BUDGET_MS)
        return 'traverse', {'seeds': [seed], 'strategy': strategy}

    out = {
        'resolve_discover_headers': ('resolve', {'sequence': ndm1, 'discover': {'kind': 'header'}}),
        'ndm1_derived_columns': traverse({'sequence': ndm1}, labels={'seed_label_kind': 'column'},
                                         bounds={'max_extension_bp': 2000}),
        'ndm1_short_radius': traverse({'sequence': ndm1}, bounds={'max_extension_bp': 10}),
    }
    for i in (0, 17):
        out[f'window_{i}'] = traverse({'sequence': records[names[i]][1000:1150]},
                                      branching={'max_label_branches': 'unlimited'},
                                      bounds={'max_extension_bp': 500})
    return out


def to_prototype(kind, request):
    """the v1 request in the prototype's schema, with the knobs v1 fixes stated"""
    if kind == 'resolve':
        return request
    st = request.get('strategy', {})
    branching = dict(st.get('branching', {}))
    branching.update({'on_reconverge': 'keep', 'hairpins': 'skip', 'max_splits_per_path': 64})
    output = dict(st.get('output', {}))
    output.update({'detail': 'full', 'continuation_bp': 0})
    strategy = {
        'schema_version': 1,
        'direction': st.get('direction', 'both'),
        'labels': st.get('labels', {}),
        'branching': branching,
        'bounds': st.get('bounds', {}),
        'output': output,
    }
    return {'seeds': request['seeds'], 'strategy': strategy}


def project_traverse(response):
    results = []
    for r in response['results']:
        if 'error' in r:
            results.append({'error': r['error'].get('cause')})
            continue
        seed = r['seed']
        out = {
            'seed': {key: seed.get(key) for key in (
                'validated_seed_id', 'seed_id_mismatch', 'length_bp', 'num_kmers', 'labels',
                'labels_from_seed', 'labels_supporting_total', 'labels_dropped')},
            'dropped_labels': [[d['label'], d['runs']] for d in seed.get('dropped_labels', [])],
            'label_dict': [[l['name'], l['kind']] for l in r['label_dict']],
            'arms': {},
        }
        for name, arm in r['arms'].items():
            segments = []
            for s in arm['segments']:
                segments.append({
                    'id': s['id'], 'parents': s['parents'], 'from_bp': s['from_bp'],
                    'length_bp': s['length_bp'], 'sequence': s.get('sequence'),
                    'labels': s['labels'], 'labels_at_end': s['labels_at_end'],
                    'label_ends': sorted([e['at_bp'], e['label'], e['reason']]
                                         for e in s['events'] if e['type'] == 'label_end'),
                })
            out['arms'][name] = {
                'status': arm['status'],
                'complete_to_bp': arm['complete_to_bp'],
                'cap': arm['cap_trigger']['reason'] if 'cap_trigger' in arm else None,
                'segments': segments,
                'splits': [[s['at_bp'], s['segment'], s['children'], s['kind']]
                           for s in arm['splits']],
                'paths': [{'segments': p['segments'], 'length_bp': p['length_bp'],
                           'path_reason': p.get('path_reason'),
                           'end_labels': sorted([e['label'], e['branches']]
                                                for e in p['end_labels'])}
                          for p in arm['paths']],
            }
        results.append(out)
    return {'results': results}


def project_resolve(response):
    out = {key: response[key] for key in ('k', 'num_kmers', 'graph_runs')}
    cut = response['labels_truncated']
    out['labels_truncated'] = cut and {'kept': cut['kept'], 'total': cut['total']}
    # the labels by support, ties by name: the ranking breaks ties its own way
    out['labels'] = sorted(([l['label'], l['kind'], l['kmers_supported'], l['runs']]
                            for l in response['labels']), key=lambda l: (-l[2], l[0]))
    return out


def project(kind, response):
    return project_resolve(response) if kind == 'resolve' else project_traverse(response)


def run_prototype(metagraph, index, kind, request):
    """the prototype's answer, through its `traverse` command"""
    with tempfile.NamedTemporaryFile('w', suffix='.json', delete=False) as f:
        json.dump(to_prototype(kind, request), f)
        path = f.name
    try:
        command = [metagraph, 'traverse'] + (['--resolve'] if kind == 'resolve' else []) + [
            '--traverse-max-time-ms', str(TIME_BUDGET_MS),
            '-i', os.path.join(index, GRAPH), '-a', os.path.join(index, ANNOTATION), path]
        res = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                             timeout=TIME_BUDGET_MS / 1000 + 60)
        if res.returncode:
            sys.exit(f'{" ".join(command)} failed:\n{res.stderr}')
        return project(kind, json.loads(res.stdout))
    finally:
        os.unlink(path)


def start_server(metagraph, index):
    """a server on the index: (process, its address)"""
    with socket.socket() as s:
        s.bind(('127.0.0.1', 0))
        port = s.getsockname()[1]
    server = subprocess.Popen(
        [metagraph, 'server_query', '-i', os.path.join(index, GRAPH),
         '-a', os.path.join(index, ANNOTATION), '--address', '127.0.0.1', '--port', str(port),
         '--traverse-max-time-ms', str(TIME_BUDGET_MS)],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    address = f'http://127.0.0.1:{port}'
    # the server answers 503 while it loads the index
    deadline = time.time() + 300
    while True:
        try:
            with urllib.request.urlopen(address + '/stats', timeout=5):
                return server, address
        except (urllib.error.URLError, OSError):
            pass
        if server.poll() is not None or time.time() > deadline:
            server.kill()
            sys.exit('The server did not get ready')
        time.sleep(0.2)


def post(address, kind, request):
    req = urllib.request.Request(f'{address}/{kind}', data=json.dumps(request).encode())
    try:
        with urllib.request.urlopen(req, timeout=TIME_BUDGET_MS / 1000 + 60) as res:
            return project(kind, json.loads(res.read()))
    except urllib.error.HTTPError as e:
        sys.exit(f'/{kind} failed with {e.code}:\n{e.read().decode()}')


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--metagraph', required=True)
    parser.add_argument('--index', required=True, help='the mini RefSeq index directory')
    parser.add_argument('--record', action='store_true',
                        help='write the golden files from this binary (the prototype)')
    parser.add_argument('--golden-dir', default=GOLDEN_DIR)
    args = parser.parse_args()

    server, address = (None, None) if args.record else start_server(args.metagraph, args.index)
    try:
        failed = compare(args, address)
    finally:
        if server:
            server.kill()
    if failed:
        sys.exit(f'{len(failed)} request(s) differ from the prototype: {", ".join(failed)}')


def compare(args, address):
    failed = []
    for name, (kind, request) in sorted(requests(args.index).items()):
        got = (run_prototype(args.metagraph, args.index, kind, request) if args.record
               else post(address, kind, request))
        golden = os.path.join(args.golden_dir, name + '.json')
        if args.record:
            os.makedirs(args.golden_dir, exist_ok=True)
            with open(golden, 'w') as f:
                json.dump(got, f, sort_keys=True, separators=(',', ':'))
                f.write('\n')
            print(f'recorded {name}')
            continue
        with open(golden) as f:
            expected = json.load(f)
        if got != expected:
            failed.append(name)
            print(f'DIFFERS: {name}')
        else:
            print(f'same: {name}')
    return failed


if __name__ == '__main__':
    main()
