#!/usr/bin/env python3
"""
Parity of `metagraph traverse` with the labeled-traversal prototype (gr/labeled-traversal,
d5d26baa) on the v1 feature set.

The requests (REQUESTS below) are v1 requests on the mini RefSeq index
(fetch_mini_refseq.sh). Both binaries' responses are projected onto the fields v1 states;
the prototype's projections are committed as golden files (tests/data/traverse_parity/).

    # check this build against the golden files (CI)
    traverse_parity.py --metagraph build/metagraph --index build/mini_refseq

    # record the golden files with a build of the prototype
    traverse_parity.py --metagraph <d5d26baa>/build/metagraph --index build/mini_refseq --record

The prototype is sent each request translated into its schema with the knobs v1 fixes stated
explicitly: no reconvergence merging (paths are kept apart), hairpins skipped, change cost
forbid, 64 splits per path, no continuation sequences, detail full.
"""
import argparse
import json
import os
import subprocess
import sys
import tempfile

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
        'resolve_select_max_support': ('resolve', {'sequence': ndm1, 'discover': {'kind': 'header'},
                                                   'select': {'policy': 'max_support', 'max_seeds': 3}}),
        'resolve_select_longest': ('resolve', {'sequence': ndm1, 'discover': {'kind': 'column'},
                                               'select': {'policy': 'longest_first', 'max_seeds': 2,
                                                          'max_labels_per_seed': 3}}),
        'ndm1_derived': traverse({'sequence': ndm1}, bounds={'max_extension_bp': 3000}),
        'ndm1_derived_columns': traverse({'sequence': ndm1}, labels={'seed_label_kind': 'column'},
                                         bounds={'max_extension_bp': 2000}),
        'ndm1_branches_2': traverse({'sequence': ndm1}, branching={'max_label_branches': 2},
                                    bounds={'max_extension_bp': 3000}),
        'ndm1_unlimited': traverse({'sequence': ndm1}, branching={'max_label_branches': 'unlimited'},
                                   bounds={'max_extension_bp': 1000}),
        'ndm1_right': traverse({'sequence': ndm1}, direction='right', bounds={'max_extension_bp': 2000}),
        'ndm1_left_no_sequences': traverse({'sequence': ndm1}, direction='left',
                                           bounds={'max_extension_bp': 2000},
                                           output={'sequences': False}),
        'ndm1_capped_labels': traverse({'sequence': ndm1}, labels={'max_seed_labels': 5},
                                       bounds={'max_extension_bp': 1000}),
        'ndm1_short_radius': traverse({'sequence': ndm1}, bounds={'max_extension_bp': 10}),
        'ndm1_prefix': traverse({'sequence': ndm1[:261]}, bounds={'max_extension_bp': 1500}),
        'ndm1_explicit_with_dropped': traverse(
            {'seed_id': 'x', 'sequence': ndm1,
             'labels': ['NZ_CM008882.1', 'NZ_MOIC01000010.1', '562']},
            bounds={'max_extension_bp': 1500}),
    }
    for i in (0, 17, 33):
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
                'labels_from_seed', 'labels_supporting_total', 'labels_dropped',
                'labels_dropped_digest')},
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
    out = {key: response[key] for key in ('k', 'num_kmers', 'graph_runs', 'labels_truncated')}
    out['labels'] = [[l['label'], l['kind'], l['kmers_supported'], l['runs']]
                     for l in response['labels']]
    out['candidates'] = [[c['kmer_interval'], c['labels']] for c in response['candidates']]
    if 'selection' in response:
        sel = response['selection']
        out['selection'] = {
            'counts': sel['counts'],
            'seeds': [[s['seed_id'], s['sequence'], s['kmer_interval'], s['labels'],
                       s['label_population']['supporting_total'],
                       s['label_population']['dropped_count'],
                       s['label_population'].get('dropped_digest'), s['overlaps_with']]
                      for s in sel['seeds']],
        }
    return out


def run(metagraph, index, kind, request, prototype):
    with tempfile.NamedTemporaryFile('w', suffix='.json', delete=False) as f:
        json.dump(to_prototype(kind, request) if prototype else request, f)
        path = f.name
    try:
        command = [metagraph, 'traverse'] + (['--resolve'] if kind == 'resolve' else []) + [
            '--traverse-max-time-ms', str(TIME_BUDGET_MS),
            '-i', os.path.join(index, GRAPH), '-a', os.path.join(index, ANNOTATION), path]
        res = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if res.returncode:
            sys.exit(f'{" ".join(command)} failed:\n{res.stderr}')
        response = json.loads(res.stdout)
    finally:
        os.unlink(path)
    return project_resolve(response) if kind == 'resolve' else project_traverse(response)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--metagraph', required=True)
    parser.add_argument('--index', required=True, help='the mini RefSeq index directory')
    parser.add_argument('--record', action='store_true',
                        help='write the golden files from this binary (the prototype)')
    parser.add_argument('--golden-dir', default=GOLDEN_DIR)
    args = parser.parse_args()

    failed = []
    for name, (kind, request) in sorted(requests(args.index).items()):
        got = run(args.metagraph, args.index, kind, request, prototype=args.record)
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
    if failed:
        sys.exit(f'{len(failed)} request(s) differ from the prototype: {", ".join(failed)}')


if __name__ == '__main__':
    main()
