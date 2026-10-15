"""One-property corruptions of an answer (ORACLE-DESIGN.md section 6).

Each mutation names the check families that must catch it (a prefix of a check name: a failing
check whose name starts with one of them counts as caught). A mutation keeps the answer
internally consistent wherever it can (`returned` follows `results`), so that it is caught by
the property it breaks and not by an incidental field.

mutations_for_entry() yields (name, families, mutated entry) for every corruption applicable to
an answered /pattern entry; incomplete=True selects the corruptions of a legitimately incomplete
answer (triplet member 3). mutations_for_resolve() does the same for a /resolve answer.
"""

import copy

from . import oracle as O


def _results(entry):
    r = entry.get('results')
    return r if isinstance(r, list) else []


def _strand_key(entry):
    rs = _results(entry)
    if rs:
        return 'strand' if 'strand' in rs[0] else 'orientation'
    return 'strand'


def _flip(o):
    return {'+': '-', '-': '+', '=': '+', 'forward': 'reverse', 'reverse': 'forward',
            'palindromic': 'forward'}.get(o, '+')


def _count_paths(entry):
    """Every count object of counts.contexts (the total and its partitions)."""
    c = entry.get('counts', {}).get('contexts')
    if not isinstance(c, dict):
        return []
    out = [c]
    if isinstance(c.get('suffix'), dict):
        out.append(c['suffix'])
    for key in ('by_offset', 'by_strand', 'by_orientation', 'by_mismatches'):
        if isinstance(c.get(key), dict):
            out += [v for v in c[key].values() if isinstance(v, dict)]
    return out


def _double(count):
    for f in ('value', 'lower', 'upper', 'estimate'):
        if isinstance(count.get(f), int):
            count[f] *= 2


def _non_context_kmer(index, exp):
    """A k-mer of the index that holds no context of the pattern (for an extra result)."""
    taken = {c.kmer for c in exp.ctxs}
    for x in index.ordered:
        if x not in taken:
            return x
    return None


def mutations_for_entry(index, config, exp, entry, incomplete=False):
    if exp is None or exp.pattern.slot or 'error' in entry or 'counts' not in entry:
        return
    rs = _results(entry)
    sk = _strand_key(entry)
    L = exp.pattern.length
    total = exp.counts['total']
    retrieval = exp.mode != 'count'

    def mutant():
        # a structured copy: the results list is copied shallowly (an answer can hold tens of
        # thousands); a mutation that changes a result copies that result first (_own)
        e = dict(entry)
        for key in ('counts', 'by_label', 'withheld', 'cut', 'stop', 'notes'):
            if key in entry:
                e[key] = copy.deepcopy(entry[key])
        if isinstance(entry.get('results'), list):
            e['results'] = list(entry['results'])
        return e

    def own(e, i):
        e['results'][i] = copy.deepcopy(e['results'][i])
        return e['results'][i]

    # soundness: a result that is no context of the pattern
    if retrieval and isinstance(entry.get('results'), list):
        x = _non_context_kmer(index, exp)
        if x is not None:
            e = mutant()
            fake = dict(rs[0]) if rs else {'node': 1, 'row': 0}
            fake.update({'kmer': x, 'offset': 0, 'instance': x[:L],
                         sk: O.STRAND_NAMES[config.mode]['fwd']})
            if exp.m:
                fake.setdefault('mismatches', 0)
                fake.setdefault('substitutions', [])
            e['results'].append(fake)
            e['returned'] = len(e['results'])
            yield 'extra_result', ('results.sound', 'results.order', 'results.node'), e
    if incomplete:
        # an incomplete answer that claims completeness, or an at_least above the truth
        if retrieval:
            e = mutant()
            e['retrieval_complete'] = True
            yield 'claim_complete', ('release',), e
        c = entry['counts']['contexts']
        if c.get('relation') in ('at_least', 'unknown'):
            e = mutant()
            e['counts']['contexts'].update({'relation': 'at_least', 'value': total + 1})
            yield 'at_least_above_truth', ('counts.total',), e
        if rs:
            e = mutant()
            e['results'].append(copy.deepcopy(rs[0]))
            e['returned'] = len(e['results'])
            yield 'duplicate_result', ('results.unique', 'results.order'), e
        return

    if rs:
        e = mutant()
        e['results'].pop(len(rs) // 2)
        e['returned'] = len(e['results'])
        yield 'drop_result', ('results.complete',), e
        e = mutant()
        e['results'].append(copy.deepcopy(rs[len(rs) // 2]))
        e['returned'] = len(e['results'])
        yield 'duplicate_result', ('results.unique', 'results.order'), e
        e = mutant()
        r = own(e, 0)
        r['offset'] = r['offset'] + 1 if r['offset'] + L < index.k else r['offset'] - 1
        yield 'offset_shift', ('results',), e
        e = mutant()
        own(e, 0)[sk] = _flip(e['results'][0][sk])
        yield 'strand_flip', ('results',), e
        e = mutant()
        inst = own(e, 0)['instance']
        e['results'][0]['instance'] = ('C' if inst[0] != 'C' else 'G') + inst[1:]
        yield 'instance_base', ('results.fields',), e
        if exp.m:
            e = mutant()
            own(e, 0)['mismatches'] = e['results'][0].get('mismatches', 0) + 1
            yield 'level_bump', ('results.fields', 'results.order'), e
            with_sub = [i for i, r in enumerate(rs) if r.get('substitutions')]
            if with_sub:
                e = mutant()
                s = own(e, with_sub[0])['substitutions'][0]
                if 'position' in s:
                    s['found'] = next(b for b in 'ACGT' if b != s['found'])
                else:
                    s['codon'] = 'TGG' if s['codon'] != 'TGG' else 'ATG'
                yield 'substitution_found', ('results.fields',), e
        keys = [(r.get(sk), r.get('kmer'), r.get('offset')) for r in rs]
        if len(set(keys)) >= 2:
            e = mutant()
            i = next(i for i in range(1, len(rs)) if keys[i] != keys[0])
            e['results'][0], e['results'][i] = e['results'][i], e['results'][0]
            yield 'swap_order', ('results.order', 'results.complete'), e
        if config.mode == 'basic':
            e = mutant()
            own(e, 0)['row'] = e['results'][0].get('row', 0) + 1
            yield 'row_shift', ('results.fields',), e
    if total > 0:
        e = mutant()
        for count in _count_paths(e):
            _double(count)
        yield 'double_counts', ('counts',), e
        e = mutant()
        c = e['counts']['contexts']
        p = next((p for p, v in c.get('by_offset', {}).items() if v.get('value')), None)
        if p is not None:
            c['by_offset'][p]['value'] += 1
            yield 'offset_count_bump', ('counts.by_offset',), e
        e = mutant()
        c = e['counts']['contexts']
        key = 'by_strand' if 'by_strand' in c else 'by_orientation'
        o = next((o for o, v in c.get(key, {}).items() if v.get('value')), None)
        if o is not None:
            c[key][o]['value'] += 1
            yield 'strand_count_bump', ('counts.by_',), e
        if exp.m:
            e = mutant()
            bm = e['counts']['contexts']['by_mismatches']
            t = next((t for t, v in bm.items() if v.get('value')), None)
            if t is not None:
                bm[t]['value'] += 1
                yield 'level_count_bump', ('counts.by_mismatches',), e
        if config.masked and entry['counts']['contexts'].get('relation') == 'exact':
            e = mutant()
            e['counts']['contexts']['relation'] = 'at_least'
            yield 'exact_to_at_least', ('counts.total',), e
    else:
        e = mutant()
        e['counts']['contexts']['value'] = 1
        yield 'zero_to_one', ('counts.total',), e
    e = mutant()
    bo = e['counts']['contexts'].get('by_offset') or {}
    if bo:
        del bo[sorted(bo)[0]]
        yield 'offset_key_drop', ('counts.keys',), e
    if retrieval:
        w = entry.get('withheld')
        if w is not None:
            e = mutant()
            e['withheld'] = None
            yield 'unwithhold', ('release',), e
        rc = entry.get('retrieval_complete')
        e = mutant()
        e['retrieval_complete'] = not rc
        yield 'flip_retrieval_complete', ('release',), e
    if exp.labels == 'all' and retrieval and rs and isinstance(rs[0].get('labels'), list):
        with_labels = [i for i, r in enumerate(rs) if r.get('labels')]
        if with_labels:
            i = with_labels[0]
            e = mutant()
            own(e, i)['labels'].pop()
            e['results'][i]['labels_total'] -= 1
            yield 'label_drop', ('labels.results',), e
            e = mutant()
            present = {l['column'] for l in own(e, i)['labels']}
            other = next((c for c in index.columns if c not in present), None)
            if other:
                lab = copy.deepcopy(e['results'][i]['labels'][0])
                lab['column'] = other
                e['results'][i]['labels'].append(lab)
                e['results'][i]['labels_total'] += 1
                yield 'label_add', ('labels.results',), e
            lab = rs[i]['labels'][0]
            if lab.get('occurrence_list'):
                e = mutant()
                o = own(e, i)['labels'][0]['occurrence_list'][0]
                a, b = map(int, o['nt_coords'].split('-'))
                o['nt_coords'] = f'{a + 1}-{b + 1}'
                yield 'occurrence_coords', ('labels.results',), e
                e = mutant()
                o = own(e, i)['labels'][0]['occurrence_list'][0]
                o['seq_id'] += 1
                yield 'occurrence_seq_id', ('labels.results',), e
        bl = entry.get('by_label')
        if isinstance(bl, list) and len(bl) >= 2:
            e = mutant()
            e['by_label'][0], e['by_label'][1] = e['by_label'][1], e['by_label'][0]
            yield 'by_label_swap', ('labels.by_label',), e
        lab = entry['counts'].get('labels') or {}
        if isinstance(lab.get('value'), int):
            e = mutant()
            e['counts']['labels']['value'] += 1
            yield 'labels_count_bump', ('labels.counts',), e


def mutations_for_resolve(answer):
    labels = answer.get('labels') if isinstance(answer, dict) else None
    if not isinstance(labels, list):
        return

    def mutant():
        return copy.deepcopy(answer)
    with_runs = [i for i, l in enumerate(labels) if l.get('runs')]
    if with_runs:
        i = with_runs[0]
        e = mutant()
        run = e['labels'][i]['runs'][-1]
        if run[1] - run[0] > 1:
            run[1] -= 1
        else:
            run[1] += 1
        e['labels'][i]['kmers_supported'] = sum(b - a for a, b in e['labels'][i]['runs'])
        yield 'run_end_shift', ('resolve.labels',), e
        e = mutant()
        dropped = e['labels'][i]['runs'].pop()
        e['labels'][i]['kmers_supported'] -= dropped[1] - dropped[0]
        yield 'run_drop', ('resolve.labels',), e
    if labels:
        e = mutant()
        e['labels'][0]['kmers_supported'] = (e['labels'][0].get('kmers_supported') or 0) + 1
        yield 'kmers_supported_bump', ('resolve.labels',), e
    if answer.get('graph_runs'):
        e = mutant()
        e['graph_runs'][0][0] += 1
        yield 'graph_runs_shift', ('resolve.graph_runs',), e
    e = mutant()
    e['num_kmers'] = (e.get('num_kmers') or 0) + 1
    yield 'num_kmers_bump', ('resolve.num_kmers',), e
    cands = answer.get('candidates')
    if isinstance(cands, list) and len(cands) >= 2 and cands[0] != cands[1]:
        e = mutant()
        e['candidates'][0], e['candidates'][1] = e['candidates'][1], e['candidates'][0]
        yield 'candidate_swap', ('resolve.candidates',), e
    if isinstance(cands, list):
        multi = [i for i, c in enumerate(cands) if len(c.get('labels', [])) >= 2]
        if multi:
            e = mutant()
            e['candidates'][multi[0]]['labels'].pop()
            yield 'candidate_label_drop', ('resolve.candidates',), e


def caught(checks, families):
    """Whether a failing check belongs to one of |families| (name prefixes)."""
    return any(c['outcome'] == 'fail' and any(c['check'].startswith(f) for f in families)
               for c in checks.items)
