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


# -- /traverse ----------------------------------------------------------------------------

_OTHER_REASON = {'dead_end': 'label_lost', 'label_lost': 'dead_end', 'max_extension_bp': 'dead_end',
                 'branch': 'label_lost', 'edge_reuse': 'rejoined_seed',
                 'rejoined_seed': 'edge_reuse', 'loss_budget': 'label_lost',
                 'trace_break': 'label_lost', 'edge_reuse_rc': 'edge_reuse'}


def mutations_for_traverse(request, answer, position, walker=None):
    """One-property corruptions of seed |position|'s result: per arm a dropped leaf, a changed
    base, a run's reason, end, branches or loss changed, a run dropped, a label end event's
    reason changed, an event dropped, a leaf's end label dropped; in annotate mode a label added
    to or dropped from a node's set and a leaf's path reason changed; the seed's dropped labels.
    A truncated arm gets the incomplete member's corruptions: claiming completeness, and a walk
    the index does not hold."""
    res = answer['results'][position]
    if 'arms' not in res:
        return
    strategy = request.get('strategy') or {}
    radius = (strategy.get('bounds') or {}).get('max_extension_bp', 5000)

    def mutant():
        a = copy.deepcopy(answer)
        return a, a['results'][position]
    for arm, aj in res['arms'].items():
        p = f'traverse.{arm}'
        segs = aj.get('segments') or []
        children = {s['id']: [c['id'] for c in segs if s['id'] in (c.get('parents') or [])]
                    for s in segs}
        merged = any(len(s.get('parents') or []) > 1 for s in segs)
        if merged:
            for si, s in enumerate(segs):
                if any(e.get('type') == 'reconverge' for e in s.get('events') or []):
                    a, r = mutant()
                    r['arms'][arm]['segments'][si]['events'] = [
                        e for e in s['events'] if e.get('type') != 'reconverge']
                    yield f'{arm}:drop_reconverge', (f'{p}.reconverge_events',), a
                    break
        if aj.get('completeness_scope') == 'united_history' and merged:
            # merge: the displays, the merge structure and the state before the first merge are
            # asserted, the united history after it on its invariants (traverse_compare.
            # _merge_invariants)
            withseq = [i for i, s in enumerate(segs) if s.get('sequence')]
            if withseq:
                a, r = mutant()
                s = r['arms'][arm]['segments'][withseq[-1]]
                b = s['sequence'][-1]
                s['sequence'] = s['sequence'][:-1] + ('A' if b != 'A' else 'C')
                yield f'{arm}:merge_segment_base', (f'{p}.sound', f'{p}.displayed_support',
                                                    f'{p}.own_routes',
                                                    f'{p}.walks_before_merge'), a
            by_id = {s['id']: s for s in segs}
            for si, s in enumerate(segs):
                ps = s.get('parents') or []
                if len(ps) < 2:
                    continue
                heads = [len(by_id.get(q, {}).get('labels_at_end') or []) for q in ps]
                # the display through a parent bringing fewer labels (SPEC 7.1; a tie would be
                # a legal swap)
                if heads[0] > heads[-1]:
                    a, r = mutant()
                    t = r['arms'][arm]['segments'][si]
                    t['parents'] = t['parents'][::-1]
                    t['labels_via_parent'] = (t.get('labels_via_parent') or [])[::-1]
                    for e in t.get('events') or []:
                        if e.get('type') == 'reconverge':
                            e['segments'] = (e.get('segments') or [])[::-1]
                    yield f'{arm}:merge_minority_first', (f'{p}.first_parent',), a
                # a label said to enter through the first parent, which did not bring it
                via = s.get('labels_via_parent') or []
                head0 = set(by_id.get(ps[0], {}).get('labels_at_end') or [])
                moved = [x for v in via[1:] for x in v if x not in head0]
                if moved:
                    a, r = mutant()
                    t = r['arms'][arm]['segments'][si]
                    v2 = t['labels_via_parent']
                    v2[0] = sorted(v2[0] + [moved[0]])
                    for v in v2[1:]:
                        if moved[0] in v:
                            v.remove(moved[0])
                    yield f'{arm}:merge_via_wrong_parent', (f'{p}.labels_via_parent',), a
                break
            # the costlier entry kept at a merge (the verifier's G13): a run the merge closed and
            # its label's cheaper continuing run trade places, labels_via_parent following
            runs_ = aj.get('runs') or []
            done = False
            for si, s in enumerate(segs):
                ps, mm = s.get('parents') or [], s.get('from_bp', 0)
                if len(ps) < 2 or done:
                    continue
                for ci, r in enumerate(runs_):
                    if r.get('end_reason') is not None or r.get('segment') not in ps or \
                            r.get('to_bp') != mm:
                        continue
                    cont = [j for j, c in enumerate(runs_) if j != ci
                            and c.get('label') == r.get('label') and c.get('segment') not in ps
                            and c.get('from_bp', 0) < mm <= c.get('to_bp', 0)
                            and float(c.get('loss', 0.0)) < float(r.get('loss', 0.0))]
                    if not cont:
                        continue
                    lab = r.get('label')
                    via0 = s.get('labels_via_parent') or []
                    came = next((j for j, v in enumerate(via0) if lab in v), None)
                    if len(via0) != len(ps) or came is None:
                        continue
                    a, res2 = mutant()
                    rr = res2['arms'][arm]['runs']
                    x, y = rr[ci], rr[cont[0]]
                    seg_y, to_y, end_y = y.get('segment'), y.get('to_bp'), y.get('end_reason')
                    # the cheaper run closed on the parent it came through, the costlier one
                    # continuing in its place, and the label entering by the costlier one's parent
                    y.update({'segment': ps[came], 'to_bp': mm, 'route_bp': 0})
                    y.pop('end_reason', None)
                    x.update({'segment': seg_y, 'to_bp': to_y, 'route_bp': 0})
                    if end_y is not None:
                        x['end_reason'] = end_y
                    via = res2['arms'][arm]['segments'][si]['labels_via_parent']
                    via[came].remove(lab)
                    k = ps.index(r['segment'])
                    via[k] = sorted(via[k] + [lab])
                    yield f'{arm}:merge_costlier_entry', (f'{p}.merge_min_entry',), a
                    done = True
                    break
            # a label_lost end after the first merge reported dead_end (the verifier's E3b): the
            # node has a successor, so it is no graph tip
            m0_ = min(s.get('from_bp', 0) for s in segs if len(s.get('parents') or []) > 1)
            late = [(si, ei) for si, s in enumerate(segs) for ei, e in enumerate(s.get('events') or [])
                    if e.get('type') == 'label_end' and e.get('reason') == 'label_lost'
                    and e.get('at_bp', 0) >= m0_]
            if late:
                a, r = mutant()
                si, ei = late[0]
                e = r['arms'][arm]['segments'][si]['events'][ei]
                e['reason'] = 'dead_end'
                for x in r['arms'][arm]['runs']:
                    if x.get('label') == e.get('label') and x.get('to_bp') == e.get('at_bp') and \
                            x.get('end_reason') == 'label_lost':
                        x['end_reason'] = 'dead_end'
                yield f'{arm}:merge_late_dead_end', (f'{p}.end_reasons_after_merge',), a
            # a run's end_reason after the first merge set to dead_end, its event unchanged (the
            # verifier's round-4 E3b and G8 miss: a superseded event under a dead_end run)
            late_runs = [i for i, x in enumerate(aj.get('runs') or [])
                         if x.get('to_bp', 0) >= m0_ and x.get('end_reason') == 'label_lost'
                         and any(e.get('type') == 'label_end' and e.get('label') == x.get('label')
                                 and e.get('at_bp') == x.get('to_bp')
                                 for s in segs if s['id'] == x.get('segment')
                                 for e in s.get('events') or [])]
            if late_runs:
                a, r = mutant()
                r['arms'][arm]['runs'][late_runs[0]]['end_reason'] = 'dead_end'
                yield f'{arm}:merge_late_run_reason', (f'{p}.run_reasons',), a
            # a run ended before the first merge, its reason changed (the keep walker's state)
            m0 = min(s.get('from_bp', 0) for s in segs if len(s.get('parents') or []) > 1)
            early = [i for i, x in enumerate(aj.get('runs') or [])
                     if x.get('to_bp', 0) < m0 and x.get('end_reason') in _OTHER_REASON]
            if early:
                a, r = mutant()
                x = r['arms'][arm]['runs'][early[0]]
                x['end_reason'] = _OTHER_REASON[x['end_reason']]
                yield f'{arm}:merge_early_reason', (f'{p}.claims_before_merge',), a
            # a run handed to a label that does not hold its end node (a label holding it too
            # would be an equivalent mutant: merge answers do not judge which co-carrier
            # carries a run)
            runs = aj.get('runs') or []
            if runs and walker is not None:
                x0 = runs[0]
                seed = request['seeds'][position]['sequence'].upper()
                walks = {}
                for s in segs:
                    seq = s.get('sequence') or ''
                    out = seq if arm == 'right' else seq[::-1]
                    par = (s.get('parents') or [None])[0]
                    walks[s['id']] = walks.get(par, '') + out
                walk = walks.get(x0.get('segment'), '')[:x0.get('to_bp', 0)]
                node = walker.walk_node(seed, arm, walk, len(walk))
                for i, d in enumerate(res.get('label_dict') or []):
                    lab = walker.label(d.get('name'))
                    if i != x0.get('label') and lab is not None and not walker.holds(lab, node):
                        a, r = mutant()
                        r['arms'][arm]['runs'][0]['label'] = i
                        yield f'{arm}:merge_run_label', (f'{p}.displayed_support',), a
                        break
            continue
        if aj.get('status') != 'complete':
            a, r = mutant()
            r['arms'][arm].update({'status': 'complete', 'complete_to_bp': radius})
            yield f'{arm}:claim_complete', (p,), a
            a, r = mutant()
            root = r['arms'][arm]['segments'][0]
            nid = max(s['id'] for s in segs) + 1
            r['arms'][arm]['segments'].append({
                'id': nid, 'parents': [root['id']],
                'from_bp': root.get('from_bp', 0) + root.get('length_bp', 0), 'length_bp': 1,
                'sequence': 'N', 'events': [], 'labels': [], 'labels_at_end': []})
            yield f'{arm}:foreign_walk', (p,), a
            continue
        leaves = [x for x in aj.get('paths') or [] if x.get('segments')]
        if leaves:
            leaf = leaves[-1]
            last = leaf['segments'][-1]
            if not children.get(last) and len(leaf['segments']) > 1:
                a, r = mutant()
                ra = r['arms'][arm]
                ra['paths'] = [x for x in ra['paths'] if x is not None and x.get('id') != leaf['id']]
                ra['segments'] = [s for s in ra['segments'] if s['id'] != last]
                ra['runs'] = [x for x in ra.get('runs') or [] if x.get('segment') != last]
                yield f'{arm}:drop_leaf', (p,), a
        withseq = [i for i, s in enumerate(segs) if s.get('sequence')]
        if withseq:
            a, r = mutant()
            s = r['arms'][arm]['segments'][withseq[-1]]
            b = s['sequence'][-1]
            s['sequence'] = s['sequence'][:-1] + ('A' if b != 'A' else 'C')
            yield f'{arm}:segment_base', (p,), a
        runs = aj.get('runs') or []
        if runs:
            i = len(runs) // 2
            a, r = mutant()
            x = r['arms'][arm]['runs'][i]
            x['end_reason'] = _OTHER_REASON.get(x['end_reason'], 'dead_end')
            yield f'{arm}:run_reason', (f'{p}.claims', f'{p}.claims_'), a
            if runs[i].get('to_bp', 0) > 0:
                a, r = mutant()
                r['arms'][arm]['runs'][i]['to_bp'] -= 1
                yield f'{arm}:run_end', (f'{p}.claims', f'{p}.shape'), a
            a, r = mutant()
            r['arms'][arm]['runs'][i]['branches'] += 1
            yield f'{arm}:run_branches', (f'{p}.claims',), a
            a, r = mutant()
            r['arms'][arm]['runs'][i]['loss'] = float(r['arms'][arm]['runs'][i]['loss']) + 0.5
            yield f'{arm}:run_loss', (f'{p}.claims',), a
            a, r = mutant()
            del r['arms'][arm]['runs'][i]
            yield f'{arm}:drop_run', (f'{p}.claims',), a
        for si, s in enumerate(segs):
            for ei, e in enumerate(s.get('events') or []):
                if e.get('type') == 'label_end':
                    a, r = mutant()
                    ev = r['arms'][arm]['segments'][si]['events'][ei]
                    ev['reason'] = _OTHER_REASON.get(ev['reason'], 'dead_end')
                    yield f'{arm}:label_end_reason', (f'{p}.label_end_events',), a
                    break
            else:
                continue
            break
        for si, s in enumerate(segs):
            idx = [ei for ei, e in enumerate(s.get('events') or [])
                   if e.get('type') in ('blocked', 'hairpin')]
            if idx:
                a, r = mutant()
                del r['arms'][arm]['segments'][si]['events'][idx[0]]
                yield f'{arm}:drop_event', (f'{p}.events',), a
                break
        with_end = [i for i, x in enumerate(aj.get('paths') or []) if x.get('end_labels')]
        if with_end:
            a, r = mutant()
            r['arms'][arm]['paths'][with_end[0]]['end_labels'].pop()
            yield f'{arm}:drop_end_label', (f'{p}.end_labels',), a
        for si, s in enumerate(segs):
            sets = s.get('label_sets') or []
            if sets:
                a, r = mutant()
                ls = r['arms'][arm]['segments'][si]['label_sets'][0]
                if ls.get('labels'):
                    ls['labels'] = ls['labels'][:-1]
                    yield f'{arm}:node_label_drop', (f'{p}.node_labels',), a
                a, r = mutant()
                ls = r['arms'][arm]['segments'][si]['label_sets'][0]
                ls['labels'] = sorted(set(ls.get('labels') or []) | {len(res.get('label_dict')
                                                                         or [])})
                yield f'{arm}:node_label_add', (f'{p}.node_labels',), a
                break
        reasons = [i for i, x in enumerate(aj.get('paths') or []) if x.get('path_reason')]
        if reasons and (strategy.get('labels') or {}).get('mode') == 'annotate':
            a, r = mutant()
            x = r['arms'][arm]['paths'][reasons[0]]
            x['path_reason'] = _OTHER_REASON.get(x['path_reason'], 'dead_end')
            yield f'{arm}:path_reason', (f'{p}.leaves',), a
    # multiplicity and order (verification H1): a listing twice, a successor enumerated twice
    for arm, aj in res['arms'].items():
        p = f'traverse.{arm}'
        segs = aj.get('segments') or []
        if aj.get('runs'):
            a, r = mutant()
            r['arms'][arm]['runs'].append(copy.deepcopy(aj['runs'][-1]))
            yield f'{arm}:duplicate_run', (f'{p}.unique', f'{p}.claims'), a
        paths = [x for x in aj.get('paths') or [] if len(x.get('segments') or []) >= 2]
        if paths:
            pth = paths[-1]
            last = next(s for s in segs if s['id'] == pth['segments'][-1])
            a, r = mutant()
            ra = r['arms'][arm]
            nid = max(s['id'] for s in segs) + 1
            clone = copy.deepcopy(last)
            clone['id'] = nid
            ra['segments'].append(clone)
            runs = ra.get('runs') or []
            remap = {}
            for i, x in enumerate(list(runs)):
                if x.get('segment') == last['id']:
                    y = copy.deepcopy(x)
                    y['segment'] = nid
                    remap[i] = len(runs)
                    runs.append(y)
            pc = copy.deepcopy(pth)
            pc['id'] = max(x['id'] for x in ra['paths']) + 1
            pc['segments'] = pth['segments'][:-1] + [nid]
            for e in pc.get('end_labels') or []:
                if e.get('run') in remap:
                    e['run'] = remap[e['run']]
            ra['paths'].append(pc)
            yield f'{arm}:duplicate_subtree', (f'{p}.unique',), a
        for si, s in enumerate(segs):
            if s.get('events'):
                a, r = mutant()
                ev = r['arms'][arm]['segments'][si]['events']
                ev.append(copy.deepcopy(ev[0]))
                ev.sort(key=lambda e: e.get('at_bp', 0))
                yield f'{arm}:duplicate_event', (f'{p}.unique',), a
                break
        for si, s in enumerate(segs):
            evs = s.get('events') or []
            pos = [e.get('at_bp', 0) for e in evs]
            if len(set(pos)) >= 2:
                a, r = mutant()
                r['arms'][arm]['segments'][si]['events'].reverse()
                yield f'{arm}:event_order', (f'{p}.order',), a
                break
        for xi, sp in enumerate(aj.get('splits') or []):
            if len(sp.get('branches') or []) >= 2:
                a, r = mutant()
                r['arms'][arm]['splits'][xi]['branches'].reverse()
                yield f'{arm}:split_order', (f'{p}.order',), a
                break
    dropped = (res.get('seed') or {}).get('dropped_labels') or []
    if dropped:
        a, r = mutant()
        r['seed']['dropped_labels'] = dropped[1:]
        yield 'seed:dropped', ('traverse.seed.dropped',), a
