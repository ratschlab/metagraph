#!/usr/bin/env python3
"""Run the oracle suite (ORACLE-DESIGN.md section 7): build the fixture indexes with the engine,
serve each on 127.0.0.1, run the generated requests, compare every answer with the oracle, apply
the mutations, write oracle_report.json and oracle_report.md, and stop every server.

    python3 metagraph/scripts/suite/oracle/run_oracle.py --metagraph PATH --out DIR
        [--seed N] [--configs basic_k11,canonical_k11,...] [--requests N] [--budget N]
        [--resolve N] [--no-mutations]

DIR must not exist or be empty (nothing is ever deleted). Exit status: 0 when every check
passed (inconclusive and exploratory outcomes allowed), 1 when one failed, 2 on a set-up error
(a build or a server that did not start).
"""

import argparse
import datetime
import gzip
import json
import os
import signal
import subprocess
import sys
import traceback

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from oracle import fixtures, indexes, report  # noqa: E402
from oracle import compare as C  # noqa: E402
from oracle import mutations as M  # noqa: E402
from oracle import oracle as O  # noqa: E402
from oracle import requests_gen  # noqa: E402
from oracle.server import Server, ServerError  # noqa: E402

DEFAULT_CONFIGS = list(indexes.CONFIGS)
SERVER_MAX_CONTEXTS = 100000

_running = []


def _stop_all(*_):
    for s in list(_running):
        s.stop()
    _running.clear()


def _signal(signum, frame):
    _stop_all()
    sys.exit(128 + signum)


def utc():
    return datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')


def config_for(name):
    mode, k, masked, anno = indexes.CONFIGS[name]
    return C.Config(name, mode, k, masked=masked,
                    annotation='budgeted' if anno == 'row_diff_brwt_coord' else 'unbudgeted')


def run_config(name, built, mg, out, columns, args, rep):
    cfg = config_for(name)
    index = O.Index(fixtures.all_records(columns), cfg.k, cfg.mode)
    reqs = requests_gen.generate(index, name, args.seed, labels_ok=True,
                                 resolve_ok=cfg.mode == 'basic',
                                 header_labels=cfg.mode == 'basic',
                                 n_pattern=args.requests, n_budget=args.budget,
                                 n_resolve=args.resolve)
    raw = os.path.join(out, 'raw', name)
    os.makedirs(raw, exist_ok=True)
    info = {'requests': len(reqs), 'entries': 0, 'outcomes': {}, 'server_log':
            os.path.join(out, 'logs', f'server_{name}.log')}
    rep['configs'][name] = info
    server = Server(mg, built['server_args'], info['server_log'])
    _running.append(server)
    try:
        status, caps = server.get('pattern/capabilities')
        if status == 200 and isinstance(caps, dict):
            info.update({'graph_mode': caps.get('graph_mode'), 'k': caps.get('k'),
                         'counting': caps.get('counting'), 'mask': caps.get('mask'),
                         'placement': caps.get('placement'), 'annotation': caps.get('annotation')})
            with open(os.path.join(raw, 'pattern_capabilities.json'), 'w') as f:
                json.dump(caps, f, indent=1, sort_keys=True)
        for item in reqs:
            status, body = server.post(item['route'], item['request'])
            with gzip.open(os.path.join(raw, item['id'] + '.json.gz'), 'wt') as f:
                json.dump({'request': item['request'], 'status': status, 'answer': body}, f,
                          sort_keys=True)
            record = {'id': item['id'], 'config': name, 'route': item['route'],
                      'class': item['class'], 'status': status, 'entries': []}
            rep['requests'].append(record)
            if item['route'] == 'resolve':
                checks, want = C.compare_resolve(index, cfg, item['request'], body, status)
                _record(rep, info, record, item, None, checks, want, body, None)
                if not args.no_mutations and checks.outcome() == C.PASS and status == 200:
                    for mname, fam, mutated in M.mutations_for_resolve(body):
                        mc, _ = C.compare_resolve(index, cfg, item['request'], mutated)
                        _mutation(rep, item, None, mname, fam, M.caught(mc, fam))
                continue
            if status != 200 or not isinstance(body, dict):
                record['outcome'] = 'fail'
                record['entries'].append({'position': None, 'outcome': 'fail',
                                          'failed': [{'check': 'http', 'expected': 200,
                                                      'got': [status, body]}]})
                _bump(info, 'fail')
                rep['disagreements'].append({'config': name, 'request_id': item['id'],
                                             'position': None, 'request': item['request'],
                                             'failed': record['entries'][0]['failed']})
                continue
            max_contexts = min(item['request'].get('max_contexts', SERVER_MAX_CONTEXTS),
                               SERVER_MAX_CONTEXTS)
            compared = C.compare_pattern(index, cfg, item['request'], body, max_contexts)
            sticky = None
            for pos, checks, exp in compared:
                entry = body['patterns'][pos] if pos is not None else None
                _record(rep, info, record, item, pos, checks, exp, entry, exp)
                if args.no_mutations or pos is None or exp is None:
                    continue
                outcome = checks.outcome()
                if outcome in (C.PASS, C.INCONCLUSIVE):
                    incomplete = outcome == C.INCONCLUSIVE
                    for mname, fam, mutated in M.mutations_for_entry(index, cfg, exp, entry,
                                                                      incomplete):
                        mc, _ = C.compare_entry(index, cfg, item['request'], pos, mutated,
                                                sticky, max_contexts, exp)
                        _mutation(rep, item, pos, mname, fam, M.caught(mc, fam))
                if sticky is None:
                    sticky = C.sticky_stop(entry.get('stop'), item['request'])
            record['outcome'] = _worst([e['outcome'] for e in record['entries']])
    finally:
        server.stop()
        if server in _running:
            _running.remove(server)


def _bump(info, outcome):
    info['outcomes'][outcome] = info['outcomes'].get(outcome, 0) + 1
    info['entries'] += 1


def _worst(outcomes):
    for o in ('fail', 'inconclusive', 'exploratory', 'pass'):
        if o in outcomes:
            return o
    return 'pass'


def _record(rep, info, record, item, pos, checks, exp, entry, pattern_exp):
    outcome = checks.outcome()
    exploratory = [c for c in checks.items if c['outcome'] == C.EXPLORATORY]
    if outcome == C.PASS and exploratory:
        outcome = 'exploratory'
    failed = checks.failed()
    record['entries'].append({
        'position': pos, 'outcome': outcome, 'checks': len(checks.items),
        'failed': failed,
        'inconclusive': [c['check'] for c in checks.items if c['outcome'] == C.INCONCLUSIVE]})
    if item['route'] == 'resolve':
        record['outcome'] = outcome
    _bump(info, outcome)
    for c in exploratory:
        rep['exploratory'].append({'config': record['config'], 'request_id': item['id'],
                                   'check': c['check'], 'note': c.get('note'),
                                   'expected': c.get('expected'), 'got': c.get('got')})
    if failed:
        d = {'config': record['config'], 'request_id': item['id'], 'position': pos,
             'request': item['request'], 'failed': failed,
             'engine': report.engine_excerpt(entry)}
        if item['route'] == 'pattern':
            d['oracle'] = report.oracle_summary(pattern_exp)
        else:
            d['oracle'] = exp
        rep['disagreements'].append(d)


def _mutation(rep, item, pos, name, families, was_caught):
    m = rep['mutations']
    m['applied'] += 1
    b = m['by_mutation'].setdefault(name, {'applied': 0, 'caught': 0})
    b['applied'] += 1
    if was_caught:
        m['caught'] += 1
        b['caught'] += 1
    else:
        m['missed'].append({'request_id': item['id'], 'position': pos, 'mutation': name,
                            'families': list(families)})


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--metagraph', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--seed', type=int, default=1)
    p.add_argument('--configs', default=','.join(DEFAULT_CONFIGS))
    p.add_argument('--requests', type=int, default=60)
    p.add_argument('--budget', type=int, default=12)
    p.add_argument('--resolve', type=int, default=30)
    p.add_argument('--no-mutations', action='store_true')
    args = p.parse_args(argv)

    mg = os.path.abspath(args.metagraph)
    out = os.path.abspath(args.out)
    if any(c.isspace() for c in out):
        print('--out must not contain whitespace (column rename rules)', file=sys.stderr)
        return 2
    if os.path.exists(out) and os.listdir(out):
        print(f'{out} is not empty: choose a new --out (nothing is ever deleted)', file=sys.stderr)
        return 2
    configs = [c for c in args.configs.split(',') if c]
    unknown = [c for c in configs if c not in indexes.CONFIGS]
    if unknown:
        print(f'unknown configs: {unknown}; known: {list(indexes.CONFIGS)}', file=sys.stderr)
        return 2
    os.makedirs(os.path.join(out, 'logs'), exist_ok=True)
    signal.signal(signal.SIGINT, _signal)
    signal.signal(signal.SIGTERM, _signal)

    version = subprocess.run([mg, '--version'], capture_output=True, text=True).stdout.strip()
    rep = {'engine': {'path': mg, 'version': version}, 'seed': args.seed,
           'started_utc': utc(), 'configs': {}, 'requests': [], 'disagreements': [],
           'exploratory': [],
           'mutations': {'applied': 0, 'caught': 0, 'missed': [], 'by_mutation': {}},
           'scope': {'covered': report.COVERED, 'not_covered': report.NOT_COVERED}}
    columns = fixtures.build_records()
    code = 0
    try:
        built = indexes.build_all(mg, out, columns, configs)
        for name in configs:
            run_config(name, built[name], mg, out, columns, args, rep)
    except (indexes.BuildError, ServerError) as e:
        rep['setup_error'] = str(e)
        print(f'set-up error: {e}', file=sys.stderr)
        code = 2
    except Exception:            # a bug of the suite: still stop the servers and report
        rep['setup_error'] = traceback.format_exc()
        print(rep['setup_error'], file=sys.stderr)
        code = 2
    finally:
        _stop_all()
    totals = {'pass': 0, 'fail': 0, 'inconclusive': 0, 'exploratory': 0}
    for info in rep['configs'].values():
        for k, v in info['outcomes'].items():
            totals[k] = totals.get(k, 0) + v
    verdict = 'FAIL' if totals['fail'] or code else 'PASS'
    rep['summary'] = dict(totals, entries=sum(totals.values()), requests=len(rep['requests']),
                          verdict=verdict + (' (set-up error)' if code == 2 else ''))
    rep['ended_utc'] = utc()
    report.write(os.path.join(out, 'oracle_report.json'), os.path.join(out, 'oracle_report.md'),
                 rep)
    print(f"{rep['summary']['verdict']}: {totals} -> {os.path.join(out, 'oracle_report.md')}")
    if code:
        return code
    return 1 if totals['fail'] else 0


if __name__ == '__main__':
    sys.exit(main())
