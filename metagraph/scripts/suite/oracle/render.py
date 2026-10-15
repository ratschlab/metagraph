"""Answers rendered from the oracle in the engine's JSON shape, for the offline tests: a valid
answer (triplet member 1), and a legitimately incomplete one (member 3: a stated max_steps stop,
at_least counts, a prefix of the list). The comparison must pass the first and judge the second
inconclusive with its invariants passing; mutations.py corrupts both (member 2).

Node ids are fake but valid: the rank of the k-mer in BOSS order, plus one (row = node - 1).
"""

from . import oracle as O


def node_and_row(index, kmer):
    """Fake but valid ids (SPEC 7.10): BASIC and CANONICAL nodes in BOSS order (rank + 1), the
    row the node of the pair's smaller node minus 1 on CANONICAL; PRIMARY stores the
    co-lexicographically smaller k-mer of each pair (ids 1..D in that order), its reverse
    complement is virtual at + D, and the row is the stored k-mer's."""
    rc = O.revcomp(kmer)
    if index.mode == 'basic':
        node = index.rank[kmer] + 1
        return node, node - 1
    if index.mode == 'canonical':
        node = index.rank[kmer] + 1
        return node, min(node, index.rank[rc] + 1) - 1
    stored = getattr(index, '_stored', None)
    if stored is None:
        reps = sorted({min(x, O.revcomp(x), key=O.colex_key) for x in index.kmers},
                      key=O.colex_key)
        stored = {x: i + 1 for i, x in enumerate(reps)}
        index._stored = stored
    rep = min(kmer, rc, key=O.colex_key)
    base = stored[rep]
    return (base if kmer == rep else base + len(stored)), base - 1


def _count(value, unit='graph_contexts', relation='exact'):
    return {'value': value if relation != 'unknown' else None, 'relation': relation, 'unit': unit}


def _subs(pattern, subs):
    if pattern.kind == 'protein':
        return [{'residue': i, 'expected': e, 'found': f, 'codon': c} for i, e, f, c in subs]
    return [{'position': i, 'expected': e, 'found': f} for i, e, f in subs]


def render_entry(index, config, request, position, max_contexts=10 ** 9, stop_after=None):
    """The entry of pattern |position|. |stop_after|: render a max_steps stop in discovery that
    found only the first |stop_after| contexts (counts at_least, the list their prefix)."""
    p = request['patterns'][position]
    kind = next(k for k in ('dna', 'iupac', 'protein') if k in p)
    code = request.get('genetic_code', 1)
    pat = O.Pattern(kind, p[kind], code)
    entry = {'id': p.get('id'), 'kind': kind}
    if pat.slot:
        entry['error'] = {'code': pat.slot, 'message': 'rendered'}
        return entry
    scope = request.get('scope', 'any_offset')
    strands = request.get('strands', 'both')
    m = request.get('max_mismatches', 0)
    mode = request.get('mode', 'all_or_count')
    labels = (request.get('output') or {}).get('labels', 'none')
    if config.mode == 'primary' and scope == 'suffix':
        entry.update({'pattern': pat.text, 'length': pat.length,
                      'error': {'code': 'scope_unsupported', 'message': 'rendered'}})
        return entry
    ctxs = O.contexts(index, pat, scope, strands, m)
    truth = O.counts(index, pat, ctxs, scope, strands, m)
    stopped = stop_after is not None and not pat.no_instance
    found = ctxs
    if stopped:
        # a stop in the work of level 0 (by_mismatches: level 0 at_least, the levels above
        # unknown): the contexts found are some of level 0's
        found = [c for c in ctxs if c.level == 0][:stop_after]
    rel = 'at_least' if stopped else 'exact'

    def part(v):
        return _count(min(v, len(found)) if stopped else v, relation=rel)
    contexts = _count(len(found), relation=rel)
    contexts['suffix'] = part(truth['suffix'])
    contexts['by_offset'] = {o: _count(sum(1 for c in found if str(c.offset) == o), relation=rel)
                             for o in truth['by_offset']}
    key = 'by_strand' if config.mode == 'basic' else 'by_orientation'
    contexts[key] = {o: _count(sum(1 for c in found
                                   if ('both' if c.orientation == '=' else c.orientation) == o),
                               relation=rel)
                     for o in truth['by_orient']}
    if m:
        contexts['by_mismatches'] = {t: _count(sum(1 for c in found if str(c.level) == t),
                                               relation=rel if not stopped else 'at_least')
                                     for t in truth['by_mismatches']}
    entry.update({
        'pattern': pat.text, 'length': pat.length, 'information_bits': pat.information_bits(),
        'anchor_information_bits': None, 'min_anchor_information_bits': None,
        'mode': mode, 'scope': scope, 'strands': O.orientations(index, pat, strands),
        'palindromic': pat.palindromic, 'absence_scope':
            'suffix_only' if scope == 'suffix' else 'any_offset',
        'determinism': 'full', 'timing': {'elapsed_ms': 0.0},
        'work': {'ranges_visited': 0, 'mask_scans': 0, 'steps': 0},
        'stop': ({'phase': 'discovery', 'reason': 'max_steps'} if stopped else None),
        'notes': [],
        'counts': {'contexts': contexts, 'labels': _count(None, 'labels', 'unknown'),
                   'occurrences': _count(None, 'placed_occurrences', 'unknown')},
    })
    if m:
        entry['stop'] = entry['stop'] and dict(entry['stop'], mismatch_level=0)
        entry.update({'max_mismatches': m, 'mismatch_information_bits':
                      pat.mismatch_information_bits(m), 'min_anchor_mismatch_information_bits':
                      None})
        entry['work']['steps_by_mismatches'] = {str(t): 0 for t in range(m + 1)}
        if stopped:
            # a stop in level 0's work: level 0 at_least, the levels above unknown
            for t in range(1, m + 1):
                contexts['by_mismatches'][str(t)] = _count(None, relation='unknown')
    if kind == 'protein':
        entry['residues'] = len(pat.residues)
        entry['genetic_code'] = code
    if config.mode != 'basic':
        entry['notes'].append('strand_unknown_canonical')
    if pat.no_instance:
        entry['notes'].append('no_stop_codon')
    if mode == 'count':
        entry['retrieval_complete'] = False
        if labels == 'all':
            entry['notes'].append('annotation_not_read')
        return entry
    released, withheld, cut = found, None, None
    if stopped:
        if mode == 'all_or_count':
            released, withheld = [], {'reason': 'discovery_budget'}
        else:
            released, cut = found[:max_contexts], {'reason': 'max_steps'}
    elif len(ctxs) > max_contexts:
        if mode == 'all_or_count':
            released, withheld = [], {'reason': 'count_above_threshold'}
        else:
            released, cut = ctxs[:max_contexts], {'reason': 'max_contexts'}
    sk = 'strand' if config.mode == 'basic' else 'orientation'
    if config.mode == 'primary':
        # PRIMARY answers in wrapper-id order (SPEC 7.9): the cut is the first by node
        def pkey(c):
            return (c.level, node_and_row(index, c.kmer)[0], c.offset,
                    O.ORIENTATION_RANK[c.orientation])
        ordered = sorted(found if stopped else ctxs, key=pkey)
        if cut is not None or withheld is not None or stopped:
            released = ordered[:len(released)] if not stopped else released
        else:
            released = ordered
    results = []
    for c in released:
        node, row = node_and_row(index, c.kmer)
        r = {'kmer': c.kmer, 'instance': c.kmer[c.offset:c.offset + pat.length],
             'offset': c.offset, sk: c.orientation, 'node': node, 'row': row}
        if m:
            r['mismatches'] = c.level
            r['substitutions'] = _subs(pat, c.subs)
        results.append(r)
    complete = withheld is None and cut is None and not stopped
    entry.update({'withheld': withheld, 'cut': cut, 'returned': len(results),
                  'results': results, 'retrieval_complete': complete})
    if labels == 'all':
        entry.update({'placement': config.placement, 'annotation': config.annotation,
                      'rows_refused': [], 'anchors_truncated': [], 'labels_cut': None,
                      'occurrences_cut': None})
        if withheld is not None:
            entry['by_label'] = None
            return entry
        view = O.label_view(index, pat, released)
        for r, want in zip(results, view['results']):
            r.update({'support': 'kmer', 'labels_status': 'complete',
                      'labels_total': len(want['labels'])})
            labs = []
            for col in want['labels']:
                lab = {'column': col, 'support': 'kmer'}
                if config.placement == 'record':
                    occ = want['occurrences'][col]
                    lab['occurrences'] = _count(len(occ), 'placed_occurrences')
                    lab['occurrence_list'] = [
                        {'seq_id': o[0], 'record': o[1], 'strand': o[2], 'nt_coords': o[3],
                         'nt_length': o[4]} for o in occ]
                labs.append(lab)
            r['labels'] = labs
        lrel = 'exact' if complete else 'at_least'
        entry['by_label'] = [
            {'graph': None, 'column': b['column'], 'contexts': _count(b['contexts'], relation=lrel),
             'contexts_suffix': _count(b['contexts_suffix'], relation=lrel),
             'occurrences': (_count(b['occurrences'], 'placed_occurrences', lrel)
                             if config.placement == 'record'
                             else _count(None, 'placed_occurrences', 'unknown'))}
            for b in view['by_label']]
        entry['counts']['labels'] = _count(view['labels'], 'labels', lrel)
        if config.placement == 'record':
            entry['counts']['occurrences'] = _count(view['occurrences'], 'placed_occurrences',
                                                    lrel)
    return entry


def render_answer(index, config, request, max_contexts=10 ** 9, stop_at=None):
    """The whole /pattern answer. |stop_at|: (position, contexts found) of a max_steps stop; the
    later patterns carry the sticky stop with unknown counts."""
    entries = []
    for i in range(len(request['patterns'])):
        if stop_at is None or i < stop_at[0]:
            entries.append(render_entry(index, config, request, i, max_contexts))
        elif i == stop_at[0]:
            entries.append(render_entry(index, config, request, i, max_contexts, stop_at[1]))
        else:
            e = render_entry(index, config, request, i, max_contexts, 0)
            if 'counts' in e and 'no_stop_codon' not in e.get('notes', []):
                _unknown(e, request)
            entries.append(e)
    block = {'k': index.k, 'graph_mode': index.mode, 'alphabet': '$ACGT',
             'strand_stated': index.mode == 'basic', 'index_ns': None, 'index_fp': None,
             'release': ''}
    if not config.masked:
        block.update({'counting': 'upper_bound',
                      'dummy_fraction': {'value': 0.99, 'interval': [0.98, 1.0],
                                         'samples': 10000, 'source': 'sampled'}})
    return {'pattern_contract_version': 1, 'mode': request.get('mode', 'all_or_count'),
            'index': block, 'patterns': entries}


def _unknown(entry, request):
    c = entry['counts']['contexts']
    for count in [c, c['suffix']] + [v for key in ('by_offset', 'by_strand', 'by_orientation',
                                                   'by_mismatches')
                                     for v in (c.get(key) or {}).values()]:
        count['relation'] = 'unknown'
        count['value'] = None
    entry['stop'] = {'phase': 'discovery', 'reason': 'max_steps'}
    if request.get('max_mismatches', 0):
        entry['stop']['mismatch_level'] = 0


def render_resolve(index, request):
    want = O.resolve(index, request['sequence'], request['labels'],
                     request.get('min_block_kmers', 1), request.get('support', 'kmer'))
    return {'k': want['k'], 'regime': index.mode, 'num_kmers': want['num_kmers'],
            'support': request.get('support', 'kmer'), 'graph_runs': want['graph_runs'],
            'labels': want['labels'], 'labels_truncated': None,
            'candidates': want['candidates'], 'release': '', 'timing': {'elapsed_ms': 0.0}}
