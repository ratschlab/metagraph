#!/usr/bin/env python3
"""bench_pattern_retrieval.py -- a repeatable before/after benchmark of /pattern's labelled retrieval and of the
engine's budget cases (CLI).

Runs `metagraph pattern --json` of two binaries (a baseline and a candidate) on the same requests and indexes,
several times each, and reports per case: the wall time, the answer's elapsed_ms, its stops and counts, the
runs refused (a 503 `deadline` exits 1), and whether the two answers are identical once their timings are
removed (`timing` objects and the top-level elapsed). The cases exercise the verification of long paths, the
occurrence cap, the low-complexity diagnostic, the sink-edge skip and the extension's clock.

CASES
  repeat (built here, cached in OUT/indexes): one record of 30,000 A, k = 3, row_diff_coord with coordinates and
    its record mapping -- the repros:
      homopolymer-path   a 1,500-base IUPAC path (A..A W A, long_search paths), occurrences, cap 1, 500 ms
      cap-volume-10      ten AAA (one context of 29,998 placed occurrences each), cap 1, 1,500 ms
  dinucleotide (built here): one record of (AC)^15,000, k = 3 -- no run of consecutive coordinates, the join's
    worst case: (AC)^750 as a path (14,000+ chains through 1,498 k-mers), cap 1, 500 ms
  sinks (built here): 100,000 records A^18 + an 11-base code over AGT + TC and one T^30 CA, k = 31, no mask --
    100,000 consecutive sink edges before the k-mer T^29 CA (the graph of the unit test
    PatternUnmasked.LongSinkRunSkippedByRankAndSelect):
      sinks-c-251ms      C, partial, max_contexts 1, forward, time_budget_ms 251 (1 ms of work time)
      sinks-c            the same without a budget
  mini_refseq (build/mini_refseq, scripts/traversal/build_mini_refseq.sh), as built (no mask), the engine's
    repros:
      atg-10000          (ATG)^10,000 as paths, mode count, max_steps 1, 251 ms (the low-complexity
                         diagnostic of a 30,000-base pattern)
      protein-10000      M^10,000 likewise
      half-n-<B>ms       GNNNC NNNT NNNA ... (40 bases, every other base N) as paths, mode count, budgets B of
                         130-180 ms, --pattern-finalize-ms 1 --pattern-min-information-bits 0 (an
                         extension of 331 anchors between clock readings). Kind 'sweep': whether a run is late
                         depends on the machine's speed and load; compare the refused runs per binary, and run
                         it with --runs 10 or more
  mini_refseq, the normal cases (with --pattern-build-mask): blaNDM-1's first 40 bases and the 51-base chimera
    as paths, NDM-1's first 14 residues and the insulin A chain as peptide paths, the NDM forward primer and the
    16S 515F primer as contexts; all_or_count and partial

USAGE
  bench_pattern_retrieval.py --baseline OLD_BIN [--candidate NEW_BIN] [--runs 5] [--out DIR] [--cases a,b]
  # the defaults: candidate <repo>/metagraph/build/metagraph, mini_refseq <repo>/metagraph/build/mini_refseq,
  # OUT bench_out/pattern_retrieval_<timestamp>; exit code 1 when a normal (mini_refseq) answer differs
  bench_pattern_retrieval.py --report OUT      # the summary table of a run again

OUTPUT (OUT)
  indexes/            the repeat, dinucleotide and sinks indexes (reused by later runs with --out pointing here)
  requests/*.json     the request bodies
  answers/            each binary's last answer per case
  results.jsonl       one row per (case, binary, run)
  summary.md          per case: median wall and elapsed of each binary, the speed-up, stops, refused runs,
                      identical or not
Stdlib only, Python >= 3.8.
"""
import argparse
import json
import os
import statistics
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_MG = os.path.abspath(os.path.join(HERE, '..', '..'))

NDM_40 = 'ATGGAATTGCCCAATATTATGCACCCGGTCGCGAAGCTGA'
CHIMERA = 'CGGCCTCCAGAGCACTTTGTCGTTTTTGGACGGAAAATCCCTAGAACCCCT'
NDM_PEP_LONG = 'MELPNIMHPVAKLS'
INSULIN_A = 'GIVEQCCTSICSLYQLENYCN'
NDM_F = 'GGTTTGGCGATCTGGTTTTC'
P515F = 'GTGYCAGCMGCCGCGGTAA'
HALF_N = 'GNNNCNNNTNNNANNNGNNNCNNNTNNNANNNGNNNCNNN'

# the flags of every case unless it names its own: exact counts on the small indexes, no information floor
DEFAULT_FLAGS = ['--pattern-build-mask', '--pattern-min-information-bits', '0']


def run(cmd, log):
    with open(log, 'a') as f:
        f.write('$ ' + ' '.join(cmd) + '\n')
        f.flush()
        subprocess.run(cmd, check=True, stdout=f, stderr=subprocess.STDOUT)


def build_index(mg, out, name, seq, k):
    """One record |seq| (header |name|) as a BASIC graph at |k|, row_diff_coord with coordinates and .seqs."""
    d = os.path.join(out, 'indexes', name)
    graph, anno = os.path.join(d, 'graph.dbg'), os.path.join(d, name + '.row_diff_coord.annodbg')
    if os.path.exists(graph) and os.path.exists(anno) and os.path.exists(os.path.join(d, name + '.seqs')):
        return graph, anno
    os.makedirs(os.path.join(d, 'rd'), exist_ok=True)
    fa = name + '.fa'
    with open(os.path.join(d, fa), 'w') as f:
        f.write('>%s\n%s\n' % (name, seq))
    log = os.path.join(d, 'build.log')
    run([mg, 'build', '--mode', 'basic', '--graph', 'succinct', '--state', 'stat', '-k', str(k),
         '-o', os.path.join(d, 'graph'), os.path.join(d, fa)], log)
    # relative names inside d: --anno-filename makes the file name the label
    cwd = os.getcwd()
    os.chdir(d)
    try:
        run([mg, 'annotate', '-p', '1', '-i', 'graph.dbg', '--anno-filename', '--coordinates', '-o', name, fa], log)
        for stage in (0, 1, 2):
            run([mg, 'transform_anno', '--anno-type', 'row_diff', '--row-diff-stage', str(stage), '--coordinates',
                 '-i', 'graph.dbg', '-o', 'rd/out', name + '.column.annodbg'], log)
        run([mg, 'transform_anno', '--anno-type', 'row_diff_coord', '-i', 'graph.dbg', '-o', name,
             os.path.join('rd', name + '.column.annodbg')], log)
        run([mg, 'annotate', '-i', 'graph.dbg', '--anno-filename', '--index-header-coords', '-o', name, fa], log)
    finally:
        os.chdir(cwd)
    return graph, anno


def build_plain(mg, out, name, records, k):
    """|records| (headers r0, r1, ...) as a BASIC graph at |k| with a column annotation, no coordinates."""
    d = os.path.join(out, 'indexes', name)
    graph, anno = os.path.join(d, 'graph.dbg'), os.path.join(d, name + '.column.annodbg')
    if os.path.exists(graph) and os.path.exists(anno):
        return graph, anno
    os.makedirs(d, exist_ok=True)
    fa = name + '.fa'
    with open(os.path.join(d, fa), 'w') as f:
        for i, seq in enumerate(records):
            f.write('>r%d\n%s\n' % (i, seq))
    log = os.path.join(d, 'build.log')
    cwd = os.getcwd()
    os.chdir(d)
    try:
        run([mg, 'build', '--mode', 'basic', '--graph', 'succinct', '-k', str(k), '-o', 'graph', fa], log)
        run([mg, 'annotate', '-i', 'graph.dbg', '--anno-filename', '-o', name, fa], log)
    finally:
        os.chdir(cwd)
    return graph, anno


def sink_records(num_sinks=100000):
    """As the unit test: A^18 + an 11-base code over AGT + TC per record, and T^30 CA after them."""
    records = []
    for i in range(num_sinks):
        code, x = '', i
        for _ in range(11):
            code += 'AGT'[x % 3]
            x //= 3
        records.append('A' * 18 + code + 'TC')
    records.append('T' * 30 + 'CA')
    return records


def cases(out, mg_for_build, mini):
    """(name, kind, graph, annotation, request, flags) of every case; kind 'repro', 'sweep' or 'normal'."""
    out_cases = []
    g, a = build_index(mg_for_build, out, 'repeat', 'A' * 30000, 3)
    homopolymer = 'A' * 1499 + 'W'
    out_cases.append(('homopolymer-path', 'repro', g, a, {
        'patterns': [{'iupac': homopolymer}], 'mode': 'partial', 'long_search': 'paths', 'strands': 'forward',
        'output': {'labels': 'all', 'occurrences': True}, 'time_budget_ms': 500, 'max_steps': 100000,
        'max_memory_mb': 64, 'max_occurrences_per_label': 1}, DEFAULT_FLAGS))
    out_cases.append(('cap-volume-10', 'repro', g, a, {
        'patterns': [{'dna': 'AAA'}] * 10, 'mode': 'partial', 'strands': 'forward',
        'output': {'labels': 'all', 'occurrences': True}, 'time_budget_ms': 1500,
        'max_occurrences_per_label': 1, 'max_memory_mb': 256}, DEFAULT_FLAGS))
    g, a = build_index(mg_for_build, out, 'dinucleotide', 'AC' * 15000, 3)
    out_cases.append(('dinucleotide-path', 'repro', g, a, {
        'patterns': [{'dna': 'AC' * 750}], 'mode': 'partial', 'long_search': 'paths', 'strands': 'forward',
        'output': {'labels': 'all', 'occurrences': True}, 'time_budget_ms': 500, 'max_steps': 100000,
        'max_memory_mb': 64, 'max_occurrences_per_label': 1}, DEFAULT_FLAGS))
    g, a = build_plain(mg_for_build, out, 'sinks', sink_records(), 31)
    sinks = {'patterns': [{'dna': 'C'}], 'mode': 'partial', 'max_contexts': 1, 'strands': 'forward'}
    no_mask = ['--pattern-min-information-bits', '0']
    out_cases.append(('sinks-c-251ms', 'repro', g, a, dict(sinks, time_budget_ms=251), no_mask))
    out_cases.append(('sinks-c', 'repro', g, a, sinks, no_mask))
    graph = os.path.join(mini, 'graph_k31.dbg')
    anno = os.path.join(mini, 'annotation.relaxed.relabeled.row_diff_brwt_coord.annodbg')
    if not os.path.exists(graph) or not os.path.exists(anno):
        print('no mini_refseq at %s: its cases are skipped' % mini, file=sys.stderr)
        return out_cases
    # the engine's repros run on the mini as built, without its mask
    for name, pattern in (('atg-10000', {'dna': 'ATG' * 10000}), ('protein-10000', {'protein': 'M' * 10000})):
        out_cases.append((name, 'repro', graph, anno, {
            'patterns': [pattern], 'long_search': 'paths', 'mode': 'count', 'max_steps': 1,
            'time_budget_ms': 251}, []))
    for budget in (130, 140, 150, 160, 170, 180):
        out_cases.append(('half-n-%dms' % budget, 'sweep', graph, anno, {
            'patterns': [{'iupac': HALF_N}], 'mode': 'count', 'long_search': 'paths', 'time_budget_ms': budget},
            ['--pattern-finalize-ms', '1', '--pattern-min-information-bits', '0']))
    normal = []
    for mode in ('all_or_count', 'partial'):
        normal.append(('ndm40-chimera-paths-' + mode, 'normal', graph, anno, {
            'patterns': [{'dna': NDM_40}, {'dna': CHIMERA}], 'mode': mode, 'long_search': 'paths',
            'output': {'labels': 'all'}}))
        normal.append(('ndm-peptide-path-' + mode, 'normal', graph, anno, {
            'patterns': [{'protein': NDM_PEP_LONG}], 'mode': mode, 'long_search': 'paths',
            'output': {'labels': 'all'}}))
        normal.append(('insulin-peptide-path-' + mode, 'normal', graph, anno, {
            'patterns': [{'protein': INSULIN_A}], 'mode': mode, 'long_search': 'paths',
            'output': {'labels': 'all'}}))
        normal.append(('primers-contexts-' + mode, 'normal', graph, anno, {
            'patterns': [{'dna': NDM_F}, {'iupac': P515F}], 'mode': mode, 'output': {'labels': 'all'}}))
    normal.append(('ndm40-paths-cap1-verified', 'normal', graph, anno, {
        'patterns': [{'dna': NDM_40}], 'mode': 'partial', 'long_search': 'paths',
        'require_support': 'record_verified', 'max_occurrences_per_label': 1, 'output': {'labels': 'all'}}))
    return out_cases + [c + (DEFAULT_FLAGS,) for c in normal]


def untimed(v):
    if isinstance(v, dict):
        return {k: untimed(x) for k, x in v.items() if k != 'timing'}
    if isinstance(v, list):
        return [untimed(x) for x in v]
    return v


def differences(x, y, path='', out=None, limit=12):
    out = [] if out is None else out
    if len(out) >= limit:
        return out
    if isinstance(x, dict) and isinstance(y, dict):
        for k in sorted(set(x) | set(y)):
            if k not in x or k not in y:
                out.append('%s/%s %s' % (path, k, 'only in candidate' if k in y else 'only in baseline'))
            else:
                differences(x[k], y[k], '%s/%s' % (path, k), out, limit)
    elif isinstance(x, list) and isinstance(y, list):
        if len(x) != len(y):
            out.append('%s: %d vs %d items' % (path, len(x), len(y)))
        for i, (a, b) in enumerate(zip(x, y)):
            differences(a, b, '%s[%d]' % (path, i), out, limit)
    elif x != y:
        out.append('%s: %s -> %s' % (path, json.dumps(x)[:80], json.dumps(y)[:80]))
    return out


def summarize(answer):
    """Per pattern: stop, retrieval_complete, occurrences count, labels count, memory_bytes."""
    if 'patterns' not in answer:
        return {'refused': answer.get('code')}
    rows = []
    for e in answer['patterns']:
        occ = e.get('counts', {}).get('occurrences', {})
        lab = e.get('counts', {}).get('labels', {})
        rows.append({'stop': e.get('stop'), 'complete': e.get('retrieval_complete'),
                     'occurrences': (occ.get('relation'), occ.get('value')),
                     'labels': (lab.get('relation'), lab.get('value')),
                     'memory_bytes': e.get('work', {}).get('memory_bytes')})
    return rows


def measure(binary, graph, anno, request_file, answer_file, flags):
    cmd = [binary, 'pattern', '--json'] + flags + ['-i', graph, '-a', anno, request_file]
    t = time.time()
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    wall = time.time() - t
    with open(answer_file, 'wb') as f:
        f.write(p.stdout)
    try:
        answer = json.loads(p.stdout)
    except ValueError:
        answer = {'code': 'non_json', 'stderr': p.stderr.decode('utf-8', 'replace')[-300:]}
    elapsed = answer.get('timing', {}).get('elapsed_ms') if isinstance(answer.get('timing'), dict) else None
    return {'rc': p.returncode, 'wall_s': wall, 'elapsed_ms': elapsed}, answer


def report(out):
    rows = [json.loads(line) for line in open(os.path.join(out, 'results.jsonl'))]
    by = {}
    for r in rows:
        by.setdefault(r['case'], {}).setdefault(r['binary'], []).append(r)
    lines = ['| case | kind | baseline wall s (median) | candidate wall s | baseline elapsed ms | candidate elapsed ms '
             '| speed-up (elapsed) | baseline stops | candidate stops | refused runs (baseline, candidate) '
             '| identical untimed |', '|' + '---|' * 11]
    for case, b in by.items():
        base, cand = b.get('baseline', []), b.get('candidate', [])
        if not base or not cand:
            continue

        def med(rs, key):
            xs = [r[key] for r in rs if r.get(key) is not None]
            return statistics.median(xs) if xs else None

        be, ce = med(base, 'elapsed_ms'), med(cand, 'elapsed_ms')

        def stops(rs):
            # the stops of every run (a sweep's differ from run to run), a refused run by its code
            seen = set()
            for r in rs:
                for e in r['summary'] if isinstance(r['summary'], list) else [r['summary']]:
                    seen.add('refused ' + str(e['refused']) if 'refused' in e else json.dumps(e.get('stop')))
            return ';'.join(sorted(seen))

        def refused(rs):
            return sum(1 for r in rs if r['rc'])

        if base[-1]['rc'] and not cand[-1]['rc']:
            speed = 'refused -> answered'
        elif cand[-1]['kind'] == 'sweep':
            # each run stops where the clock happens to be: no speed comparison
            speed = 'n/a (sweep)'
        elif stops(base) != stops(cand):
            # not the same work: a stop that fired (or no longer fires) is no speed comparison
            speed = 'stops differ'
        else:
            speed = '%.2fx' % (be / ce) if be and ce else 'n/a'
        bw, cw = med(base, 'wall_s'), med(cand, 'wall_s')
        if cand[-1]['kind'] == 'sweep':
            identical = 'n/a (sweep)'
        else:
            identical = 'yes' if cand[-1]['identical'] else 'NO: ' + '; '.join(cand[-1]['differences'][:4])
        lines.append('| %s | %s | %.3f | %.3f | %s | %s | %s | %s | %s | %d of %d, %d of %d | %s |' % (
            case, cand[-1]['kind'], bw, cw, '%.2f' % be if be is not None else 'rc %d' % base[-1]['rc'],
            '%.2f' % ce if ce is not None else 'rc %d' % cand[-1]['rc'], speed, stops(base), stops(cand),
            refused(base), len(base), refused(cand), len(cand), identical))
    text = '\n'.join(lines) + '\n'
    with open(os.path.join(out, 'summary.md'), 'w') as f:
        f.write(text)
    print(text)
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--baseline')
    ap.add_argument('--candidate', default=os.path.join(REPO_MG, 'build', 'metagraph'))
    ap.add_argument('--mini', default=os.path.join(REPO_MG, 'build', 'mini_refseq'))
    ap.add_argument('--runs', type=int, default=5)
    ap.add_argument('--out')
    ap.add_argument('--cases', help='comma-separated case names (default: all)')
    ap.add_argument('--report', metavar='OUT')
    args = ap.parse_args()
    if args.report:
        report(args.report)
        return 0
    if not args.baseline:
        ap.error('--baseline is required')
    out = os.path.abspath(args.out or os.path.join('bench_out', time.strftime('pattern_retrieval_%Y%m%d_%H%M%S')))
    for sub in ('requests', 'answers'):
        os.makedirs(os.path.join(out, sub), exist_ok=True)
    with open(os.path.join(out, 'meta.json'), 'w') as f:
        json.dump({'baseline': args.baseline, 'candidate': args.candidate, 'runs': args.runs,
                   'start': time.strftime('%Y-%m-%dT%H:%M:%S')}, f, indent=1)
    selected = set(args.cases.split(',')) if args.cases else None
    differ = False
    with open(os.path.join(out, 'results.jsonl'), 'a') as results:
        for name, kind, graph, anno, request, flags in cases(out, args.candidate, args.mini):
            if selected and name not in selected:
                continue
            request_file = os.path.join(out, 'requests', name + '.json')
            with open(request_file, 'w') as f:
                json.dump(request, f)
            last = {}
            for binary_name, binary in (('baseline', args.baseline), ('candidate', args.candidate)):
                for i in range(args.runs):
                    m, answer = measure(binary, graph, anno, request_file,
                                        os.path.join(out, 'answers', '%s.%s.json' % (name, binary_name)), flags)
                    last[binary_name] = answer
                    row = dict(m, case=name, kind=kind, binary=binary_name, run=i, summary=summarize(answer))
                    if binary_name == 'candidate' and i == args.runs - 1:
                        diff = differences(untimed(last['baseline']), untimed(answer))
                        row['identical'] = not diff
                        row['differences'] = diff
                        differ |= bool(diff) and kind == 'normal'
                    results.write(json.dumps(row) + '\n')
                    results.flush()
            print('%-34s done' % name, file=sys.stderr)
    report(out)
    return 1 if differ else 0


if __name__ == '__main__':
    sys.exit(main())
