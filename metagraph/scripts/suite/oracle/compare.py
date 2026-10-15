"""Hold an engine answer against the oracle, check by check (ORACLE-DESIGN.md section 4).

Every check has a name (its family before the first dot: shape, notes, counts, results, release,
labels, stop, resolve) and an outcome: pass, fail, inconclusive (the property cannot be judged
from this answer: a stated budget stop, a withheld list) or exploratory (a reading of the SPEC
that phase 1 states rather than enforces). A failure carries what was expected and what came.

The comparison requires complete sets and exact counts wherever the answer is complete, by its
own claim or by the SPEC's rule for its mode and budget; a legitimately incomplete answer (a
stated stop) is judged on its observable invariants only.
"""

from . import oracle as O

PASS, FAIL, INCONCLUSIVE, EXPLORATORY = 'pass', 'fail', 'inconclusive', 'exploratory'
STOP_PHASES = {'discovery', 'mask_scan', 'extraction'}
STOP_REASONS = {'max_steps', 'time', 'max_contexts', 'max_anchors'}
BUDGET_REASONS = {'max_steps', 'time'}


class Config:
    """What the comparison needs to know of the served index."""

    def __init__(self, name, mode, k, masked=True, labels=True, placement=None,
                 annotation='budgeted'):
        self.name = name
        self.mode = mode
        self.k = k
        self.masked = masked
        self.labels = labels
        self.placement = placement or ('record' if mode == 'basic' else 'none_canonical')
        self.annotation = annotation

    @property
    def strand_key(self):
        return 'strand' if self.mode == 'basic' else 'orientation'

    @property
    def by_key(self):
        return 'by_strand' if self.mode == 'basic' else 'by_orientation'


class Checks:
    def __init__(self):
        self.items = []

    def add(self, name, ok, expected=None, got=None, note=None):
        item = {'check': name, 'outcome': PASS if ok else FAIL}
        if not ok:
            item['expected'] = _clip(expected)
            item['got'] = _clip(got)
        if note:
            item['note'] = note
        self.items.append(item)
        return ok

    def mark(self, name, outcome, note=None, expected=None, got=None):
        item = {'check': name, 'outcome': outcome}
        if note:
            item['note'] = note
        if expected is not None or got is not None:
            item['expected'] = _clip(expected)
            item['got'] = _clip(got)
        self.items.append(item)

    def failed(self):
        return [c for c in self.items if c['outcome'] == FAIL]

    def outcome(self):
        outs = {c['outcome'] for c in self.items}
        if FAIL in outs:
            return FAIL
        if INCONCLUSIVE in outs:
            return INCONCLUSIVE
        return PASS


def _clip(value, limit=4000):
    text = repr(value)
    if len(text) <= limit:
        return value
    return text[:limit] + '...'


# -- counts -------------------------------------------------------------------------------

def judge_count(checks, name, count, truth, exact_required, unit='graph_contexts', upper=None):
    """One count against the true number; |upper|: U, the upper bound a bounds count of a graph
    without its mask states (SPEC section 7.4), when the oracle knows it."""
    if not isinstance(count, dict):
        return checks.add(name, False, f'a count of {truth}', count)
    rel = count.get('relation')
    ok_unit = count.get('unit') == unit
    if rel == 'exact':
        return checks.add(name, ok_unit and count.get('value') == truth
                          and 'lower' not in count and 'estimate' not in count,
                          {'value': truth, 'relation': 'exact', 'unit': unit}, count)
    if exact_required:
        return checks.add(name, False, {'value': truth, 'relation': 'exact', 'unit': unit}, count)
    if rel == 'at_least':
        v = count.get('value')
        return checks.add(name, ok_unit and isinstance(v, int) and 0 <= v <= truth,
                          f'at_least <= {truth}', count)
    if rel == 'bounds':
        lo, hi, v = count.get('lower'), count.get('upper'), count.get('value')
        ok = (ok_unit and isinstance(lo, int) and isinstance(hi, int) and v == lo
              and lo <= truth <= hi)
        if 'estimate' in count:
            ok = ok and lo <= count['estimate'] <= hi
        if upper is not None:
            ok = ok and hi == upper
        return checks.add(name, ok, f'bounds holding {truth}'
                          + (f', upper {upper}' if upper is not None else ''), count)
    if rel == 'unknown':
        return checks.add(name, ok_unit and count.get('value') is None, 'unknown: value null',
                          count)
    return checks.add(name, False, f'a relation of the vocabulary for {truth}', count)


# -- one /pattern entry -------------------------------------------------------------------

class Expectation:
    """The oracle's complete answer to one pattern of a request."""

    def __init__(self, index, request, req_pattern):
        self.kind = next(k for k in ('dna', 'iupac', 'protein') if k in req_pattern)
        self.code = request.get('genetic_code', 1)
        self.pattern = O.Pattern(self.kind, req_pattern[self.kind], self.code)
        self.scope = request.get('scope', 'any_offset')
        self.strands = request.get('strands', 'both')
        self.m = request.get('max_mismatches', 0)
        self.mode = request.get('mode', 'all_or_count')
        self.labels = (request.get('output') or {}).get('labels', 'none')
        if self.pattern.slot:
            return
        self.ctxs = O.contexts(index, self.pattern, self.scope, self.strands, self.m)
        self.counts = O.counts(index, self.pattern, self.ctxs, self.scope, self.strands, self.m)
        self.by_identity = {c.identity: c for c in self.ctxs}
        self.orientations = O.orientations(index, self.pattern, self.strands)
        self.upper = None

    def with_upper(self, index):
        """The source-dummy contexts per level (a graph without its mask): U of a count is its
        masked value plus the dummies of the levels the answer leaves in bounds (SPEC 7.4; with
        substitutions the check of few unchecked candidates makes levels exact one by one,
        21.6, and the partitions are totals over the levels)."""
        if self.upper is None and not self.pattern.slot:
            dummies = O.dummy_contexts(index, self.pattern, self.scope, self.strands, self.m)
            self.upper = {}
            for t in range(self.m + 1):
                at = [d for d in dummies if d.level == t]
                self.upper[str(t)] = O.counts(index, self.pattern, at, self.scope, self.strands,
                                              0)
        return self.upper


def _result_identity(config, r):
    return (r.get(config.strand_key), r.get('kmer'), r.get('offset'))


def _subs_json(pattern, subs):
    if pattern.kind == 'protein':
        return [{'residue': i, 'expected': e, 'found': f, 'codon': c} for i, e, f, c in subs]
    return [{'position': i, 'expected': e, 'found': f} for i, e, f in subs]


def compare_entry(index, config, request, position, entry, previous_stop=None,
                  max_contexts=None, exp=None):
    """Checks of the entry at |position| of |request|'s answer. |previous_stop|: the stop a later
    pattern must carry after an earlier pattern's budget stop (sticky, SPEC section 7.6), or
    None. |max_contexts|: the request's effective max_contexts (the server's cap applied).
    |exp|: the oracle's Expectation of this pattern when already computed (the mutations reuse
    it)."""
    checks = Checks()
    req_pattern = request['patterns'][position]
    if exp is None:
        exp = Expectation(index, request, req_pattern)
        exp.max_contexts = (max_contexts if max_contexts is not None
                            else request.get('max_contexts', 10 ** 9))
    pat = exp.pattern
    checks.add('shape.id', entry.get('id') == req_pattern.get('id'), req_pattern.get('id'),
               entry.get('id'))
    checks.add('shape.kind', entry.get('kind') == exp.kind, exp.kind, entry.get('kind'))

    if pat.slot:
        err = entry.get('error') or {}
        checks.add('slot', err.get('code') == pat.slot and 'counts' not in entry, pat.slot, err)
        return checks, exp
    if config.mode == 'primary' and exp.scope == 'suffix':
        err = entry.get('error') or {}
        checks.add('slot', err.get('code') == 'scope_unsupported' and 'counts' not in entry,
                   'scope_unsupported', err)
        return checks, exp
    if 'error' in entry:
        checks.add('slot', False, 'an answered entry', entry.get('error'))
        return checks, exp

    _shape(checks, config, exp, entry)
    stop = entry.get('stop')
    _notes(checks, config, exp, entry, stop)

    if previous_stop is not None and not pat.no_instance:
        checks.add('stop.sticky', stop == previous_stop, previous_stop, stop)
    stopped = stop is not None
    if stopped:
        checks.mark('completeness', INCONCLUSIVE,
                    note='a stated stop: completeness and exact counts not judged')
        ok = (isinstance(stop, dict) and stop.get('phase') in STOP_PHASES
              and stop.get('reason') in STOP_REASONS)
        checks.add('stop.vocabulary', ok, 'a known phase and reason', stop)
        if exp.m:
            lvl = stop.get('mismatch_level') if isinstance(stop, dict) else None
            if stop.get('phase') in ('discovery', 'mask_scan'):
                checks.add('stop.mismatch_level', isinstance(lvl, int) and 0 <= lvl <= exp.m,
                           f'0..{exp.m}', lvl)
    if pat.no_instance:
        checks.add('stop.no_instance', stop is None, None, stop)
        checks.add('notes.no_stop_codon', 'no_stop_codon' in (entry.get('notes') or []),
                   'no_stop_codon', entry.get('notes'))

    released = _release(checks, config, exp, entry, stopped)
    _counts(checks, index, config, exp, entry, stopped, released)
    results = entry.get('results')
    if exp.mode != 'count' and isinstance(results, list):
        _results(checks, config, exp, entry, results, released, stopped)
        if exp.labels == 'all':
            _labels(checks, index, config, exp, entry, results, released, stopped)
    elif exp.mode != 'count':
        checks.add('release.results', False, 'a list of results', results)
    if exp.mode == 'count' and exp.labels == 'all':
        checks.add('labels.not_read', 'annotation_not_read' in (entry.get('notes') or [])
                   and entry['counts']['labels'].get('relation') == 'unknown',
                   'annotation_not_read and unknown labels', entry.get('notes'))
    return checks, exp


def _close(a, b):
    return isinstance(b, (int, float)) and abs(a - b) <= 1e-6 * max(1.0, abs(a))


def _shape(checks, config, exp, entry):
    pat = exp.pattern
    checks.add('shape.pattern', entry.get('pattern') == pat.text, pat.text, entry.get('pattern'))
    checks.add('shape.length', entry.get('length') == pat.length, pat.length, entry.get('length'))
    checks.add('shape.palindromic', entry.get('palindromic') == pat.palindromic,
               pat.palindromic, entry.get('palindromic'))
    checks.add('shape.strands', entry.get('strands') == exp.orientations, exp.orientations,
               entry.get('strands'))
    checks.add('shape.scope', entry.get('scope') == exp.scope, exp.scope, entry.get('scope'))
    absence = 'suffix_only' if exp.scope == 'suffix' else 'any_offset'
    checks.add('shape.absence_scope', entry.get('absence_scope') == absence, absence,
               entry.get('absence_scope'))
    checks.add('shape.mode', entry.get('mode') == exp.mode, exp.mode, entry.get('mode'))
    checks.add('shape.anchor_bits', entry.get('anchor_information_bits') is None
               and entry.get('min_anchor_information_bits') is None, None,
               (entry.get('anchor_information_bits'), entry.get('min_anchor_information_bits')))
    if exp.kind == 'protein':
        checks.add('shape.residues', entry.get('residues') == len(pat.residues),
                   len(pat.residues), entry.get('residues'))
        checks.add('shape.genetic_code', entry.get('genetic_code') == exp.code, exp.code,
                   entry.get('genetic_code'))
    bits = pat.information_bits()
    checks.add('shape.information_bits', _close(bits, entry.get('information_bits')), bits,
               entry.get('information_bits'))
    if exp.m:
        checks.add('shape.max_mismatches', entry.get('max_mismatches') == exp.m, exp.m,
                   entry.get('max_mismatches'))
        mbits = pat.mismatch_information_bits(exp.m)
        checks.add('shape.mismatch_information_bits',
                   _close(mbits, entry.get('mismatch_information_bits')), mbits,
                   entry.get('mismatch_information_bits'))
    else:
        checks.add('shape.no_mismatch_fields', 'max_mismatches' not in entry
                   and 'mismatch_information_bits' not in entry, None,
                   entry.get('max_mismatches'))


def _notes(checks, config, exp, entry, stop):
    notes = entry.get('notes') or []
    if config.mode != 'basic':
        checks.add('notes.strand_unknown_canonical', 'strand_unknown_canonical' in notes,
                   'strand_unknown_canonical', notes)
    checks.add('notes.no_long_pattern_note', 'paths_later_increment' not in notes, None, notes)
    if not exp.pattern.no_instance:
        checks.add('notes.no_stop_codon_absent', 'no_stop_codon' not in notes, None, notes)
    has_estimate = 'estimate_sampled_dummy_fraction' in notes
    if config.masked:
        checks.add('notes.no_estimate_on_masked', not has_estimate, None, notes)


def _release(checks, config, exp, entry, stopped):
    """Mode fields; -> the list of contexts the answer must hold completely and in order (the
    oracle's), or None when it is not judged complete (count mode, withheld, a stop)."""
    mode = exp.mode
    total = exp.counts['total']
    max_contexts = exp.max_contexts
    rc = entry.get('retrieval_complete')
    if mode == 'count':
        absent = all(f not in entry for f in ('results', 'returned', 'withheld', 'cut'))
        checks.add('release.count_mode', absent and rc is False,
                   'no results/returned/withheld/cut, retrieval_complete false',
                   {f: entry.get(f) for f in ('returned', 'withheld', 'cut', 'retrieval_complete')})
        return None
    withheld, cut, returned = entry.get('withheld'), entry.get('cut'), entry.get('returned')
    results = entry.get('results') if isinstance(entry.get('results'), list) else []
    checks.add('release.returned', returned == len(results), len(results), returned)
    if exp.pattern.no_instance:
        checks.add('release.no_instance', withheld is None and cut is None and rc is True
                   and results == [], 'an empty complete answer',
                   {'withheld': withheld, 'cut': cut, 'retrieval_complete': rc})
        return []
    if stopped:
        checks.add('release.stopped_not_complete', rc is False, False, rc)
        if mode == 'all_or_count':
            ok = (isinstance(withheld, dict) and withheld.get('reason') in
                  ('discovery_budget', 'threshold_crossed', 'deadline') and results == [])
            checks.add('release.stopped_withheld', ok, 'withheld for the stop, nothing released',
                       withheld)
        else:
            ok = isinstance(cut, dict) and cut.get('reason') and len(results) <= max_contexts
            checks.add('release.stopped_cut', bool(ok), f'a cut, at most {max_contexts} results',
                       {'cut': cut, 'returned': len(results)})
        return None
    unmasked_bounds = (not config.masked and entry.get('counts', {}).get('contexts', {})
                       .get('relation') == 'bounds')
    if mode == 'all_or_count':
        if total <= max_contexts and not unmasked_bounds:
            checks.add('release.all_or_count', withheld is None and cut is None and rc is True,
                       'released, complete',
                       {'withheld': withheld, 'cut': cut, 'retrieval_complete': rc})
            return exp.ctxs
        if total > max_contexts:
            ok = (withheld == {'reason': 'count_above_threshold'} and results == []
                  and rc is False)
            checks.add('release.all_or_count', ok, 'withheld count_above_threshold',
                       {'withheld': withheld, 'returned': returned, 'retrieval_complete': rc})
            return None
        # a graph without its mask: U above max_contexts may withhold a count that fits
        ok = (withheld == {'reason': 'count_above_threshold'} and results == []
              and 'threshold_upper_bound' in (entry.get('notes') or []))
        checks.add('release.all_or_count', ok, 'withheld on U with threshold_upper_bound',
                   {'withheld': withheld, 'notes': entry.get('notes')})
        return None
    # partial
    if total <= max_contexts:
        checks.add('release.partial', cut is None and withheld is None and rc is True,
                   'released, complete',
                   {'withheld': withheld, 'cut': cut, 'retrieval_complete': rc})
        return exp.ctxs
    ok = cut == {'reason': 'max_contexts'} and rc is False and len(results) == max_contexts
    checks.add('release.partial', ok, f'cut max_contexts at {max_contexts}',
               {'cut': cut, 'returned': returned, 'retrieval_complete': rc})
    return exp.ctxs[:max_contexts]


def _counts(checks, index, config, exp, entry, stopped, released):
    counts = entry.get('counts')
    if not isinstance(counts, dict) or not isinstance(counts.get('contexts'), dict):
        checks.add('counts.present', False, 'counts.contexts', counts)
        return
    c = counts['contexts']
    truth = exp.counts
    up = exp.with_upper(index) if not config.masked else None
    if exp.m:
        in_bounds = [t for t, v in (c.get('by_mismatches') or {}).items()
                     if isinstance(v, dict) and v.get('relation') == 'bounds']
    else:
        in_bounds = ['0'] if c.get('relation') == 'bounds' else []

    def u(path, key=None, levels=None):
        """U of a partition: the truth plus the dummies of the levels in bounds."""
        if up is None:
            return None
        base = truth[path] if key is None else truth[path][key]
        extra = 0
        for t in (in_bounds if levels is None else levels):
            if t in up:
                extra += up[t][path] if key is None else up[t][path][key]
        return base + extra
    # on a masked graph a completed discovery is exact everywhere; without the mask bounds are
    # allowed unless the release enumerated every candidate (SPEC section 7.4)
    complete_release = released is not None and len(released) == len(exp.ctxs)
    exact_required = not stopped and (config.masked or complete_release
                                      or exp.pattern.no_instance)
    judge_count(checks, 'counts.total', c, truth['total'], exact_required, upper=u('total'))
    judge_count(checks, 'counts.suffix', c.get('suffix'), truth['suffix'], exact_required,
                upper=u('suffix'))
    bo = c.get('by_offset') or {}
    checks.add('counts.keys.by_offset', set(bo) == set(truth['by_offset']),
               sorted(truth['by_offset']), sorted(bo))
    for p, v in truth['by_offset'].items():
        if p in bo:
            judge_count(checks, f'counts.by_offset[{p}]', bo[p], v, exact_required,
                        upper=u('by_offset', p))
    key = config.by_key
    bs = c.get(key) or {}
    checks.add('counts.keys.' + key, set(bs) == set(truth['by_orient']), sorted(truth['by_orient']),
               sorted(bs))
    for o, v in truth['by_orient'].items():
        if o in bs:
            judge_count(checks, f'counts.{key}[{o}]', bs[o], v, exact_required,
                        upper=u('by_orient', o))
    other = 'by_orientation' if key == 'by_strand' else 'by_strand'
    checks.add('counts.keys.no_' + other, other not in c, None, c.get(other))
    if exp.m:
        bm = c.get('by_mismatches') or {}
        checks.add('counts.keys.by_mismatches', set(bm) == set(truth['by_mismatches']),
                   sorted(truth['by_mismatches']), sorted(bm))
        stop = entry.get('stop') or {}
        stop_level = stop.get('mismatch_level')
        for t, v in truth['by_mismatches'].items():
            if t not in bm:
                continue
            # a stop in the work of level t leaves the levels below t complete (SPEC 21.6)
            below = stop_level is not None and int(t) < stop_level
            judge_count(checks, f'counts.by_mismatches[{t}]', bm[t], v,
                        exact_required or (below and config.masked),
                        upper=None if up is None else v + up[t]['total'])
            if stop_level is not None and int(t) > stop_level:
                checks.add(f'counts.by_mismatches[{t}].unknown_above_stop',
                           bm[t].get('relation') == 'unknown', 'unknown', bm[t])
    else:
        checks.add('counts.keys.no_by_mismatches', 'by_mismatches' not in c, None,
                   c.get('by_mismatches'))


def _results(checks, config, exp, entry, results, released, stopped):
    pat = exp.pattern
    L = pat.length
    seen, extra, dup, bad_fields = set(), [], [], []
    nodes = {}
    node_conflicts = []
    for r in results:
        ident = _result_identity(config, r)
        ctx = exp.by_identity.get(ident)
        if ctx is None:
            extra.append(ident)
            continue
        if ident in seen:
            dup.append(ident)
        seen.add(ident)
        kmer, off = r.get('kmer'), r.get('offset')
        want = {'instance': kmer[off:off + L]}
        got = {'instance': r.get('instance')}
        if exp.m:
            want['mismatches'] = ctx.level
            want['substitutions'] = _subs_json(pat, ctx.subs)
            got['mismatches'] = r.get('mismatches')
            got['substitutions'] = r.get('substitutions')
        elif 'mismatches' in r or 'substitutions' in r:
            got['mismatches'] = r.get('mismatches')
        if config.mode == 'basic':
            want['row'] = r.get('node', 0) - 1 if isinstance(r.get('node'), int) else 'node - 1'
            got['row'] = r.get('row')
        if want != got:
            bad_fields.append({'identity': ident, 'expected': want, 'got': got})
        node = r.get('node')
        if kmer in nodes and nodes[kmer] != node:
            node_conflicts.append((kmer, nodes[kmer], node))
        nodes[kmer] = node
    checks.add('results.sound', not extra, 'every result an oracle context', extra[:20])
    checks.add('results.unique', not dup, 'no identity twice', dup[:20])
    checks.add('results.fields', not bad_fields, 'instance, level, substitutions, row',
               bad_fields[:10])
    inverse = {}
    for kmer, node in nodes.items():
        inverse.setdefault(node, set()).add(kmer)
    collisions = {n: sorted(v) for n, v in inverse.items() if len(v) > 1}
    checks.add('results.node', not node_conflicts and not collisions,
               'one node per k-mer, one k-mer per node',
               {'conflicts': node_conflicts[:10], 'collisions': list(collisions.items())[:10]})
    # the order (SPEC 7.9, 21.7)
    if config.mode in ('basic', 'canonical'):
        valid = [r for r in results if _result_identity(config, r) in exp.by_identity]
        keys = [_order_key(exp, config, r) for r in valid]
        checks.add('results.order', keys == sorted(keys),
                   '(level, BOSS order, offset, orientation)', _first_disorder(keys, valid, config))
        by_rank = sorted(nodes, key=O.colex_key)
        ok = [nodes[x] for x in by_rank] == sorted(nodes[x] for x in by_rank)
        checks.add('results.node_order', ok, 'node ids increasing in BOSS order of the k-mers',
                   None if ok else 'node ids out of BOSS order')
    else:
        keys = [(r.get('mismatches', 0), r.get('node'), r.get('offset'),
                 O.ORIENTATION_RANK.get(r.get(config.strand_key), 9)) for r in results]
        checks.add('results.order', keys == sorted(keys), '(level, node, offset, orientation)',
                   _first_disorder(keys, results, config))
    if released is None:
        return
    want = [c.identity for c in released]
    got = [_result_identity(config, r) for r in results]
    if config.mode == 'primary' and len(released) < len(exp.ctxs):
        # a partial cut on PRIMARY: the first max_contexts by the answer's own node ids
        ok = set(got) <= set(exp.by_identity) and len(got) == len(want)
        checks.add('results.complete', ok, f'{len(want)} contexts of the oracle set', got[:20])
        return
    got_set = set(got)
    missing = [i for i in want if i not in got_set]
    if config.mode == 'primary':
        ok = set(got) == set(want) and len(got) == len(want)
    else:
        ok = got == want
    checks.add('results.complete', ok, {'count': len(want), 'missing': missing[:20]},
               {'count': len(got), 'first_difference': _first_diff(want, got)})


def _order_key(exp, config, r):
    return (r.get('mismatches', 0) if exp.m else 0, O.colex_key(r.get('kmer')), r.get('offset'),
            O.ORIENTATION_RANK.get(r.get(config.strand_key), 9))


def _first_disorder(keys, results, config):
    for i in range(1, len(keys)):
        if keys[i] < keys[i - 1]:
            return {'at': i, 'previous': _result_identity(config, results[i - 1]),
                    'this': _result_identity(config, results[i])}
    return None


def _first_diff(want, got):
    for i, (a, b) in enumerate(zip(want, got)):
        if a != b:
            return {'at': i, 'expected': a, 'got': b}
    if len(want) != len(got):
        return {'at': min(len(want), len(got)), 'lengths': [len(want), len(got)]}
    return None


# -- labels -------------------------------------------------------------------------------

def _labels(checks, index, config, exp, entry, results, released, stopped):
    counts = entry.get('counts') or {}
    checks.add('labels.placement', entry.get('placement') == config.placement, config.placement,
               entry.get('placement'))
    checks.add('labels.annotation', entry.get('annotation') == config.annotation,
               config.annotation, entry.get('annotation'))
    if entry.get('withheld') is not None:
        checks.add('labels.withheld', entry.get('by_label') is None
                   and counts.get('labels', {}).get('relation') == 'unknown',
                   'by_label null, labels unknown',
                   {'by_label': entry.get('by_label'), 'labels': counts.get('labels')})
        return
    valid = []
    for r in results:
        ctx = exp.by_identity.get(_result_identity(config, r))
        if ctx is not None:
            valid.append((r, ctx))
    view = O.label_view(index, exp.pattern, [ctx for _, ctx in valid])
    bad = []
    for (r, ctx), want in zip(valid, view['results']):
        got_labels = r.get('labels')
        names = [l.get('column') for l in got_labels] if isinstance(got_labels, list) else None
        problem = {}
        if r.get('labels_status') != 'complete':
            problem['labels_status'] = r.get('labels_status')
        if r.get('labels_total') != len(want['labels']):
            problem['labels_total'] = (len(want['labels']), r.get('labels_total'))
        if names != want['labels']:
            problem['labels'] = (want['labels'], names)
        elif config.placement == 'record':
            for lab in got_labels:
                occ = want['occurrences'][lab['column']]
                got_occ = [(o.get('seq_id'), o.get('record'), o.get('strand'), o.get('nt_coords'),
                            o.get('nt_length')) for o in (lab.get('occurrence_list') or [])]
                cnt = lab.get('occurrences') or {}
                if (got_occ != occ or cnt.get('value') != len(occ)
                        or cnt.get('relation') != 'exact' or lab.get('support') != 'kmer'):
                    problem.setdefault('occurrences', []).append(
                        {'column': lab['column'], 'expected': occ, 'got': got_occ,
                         'count': cnt})
        else:
            # none_canonical: labels only, no occurrence fields (SPEC 14.5)
            for lab in got_labels:
                if (lab.get('support') != 'kmer' or 'occurrence_list' in lab
                        or 'occurrences' in lab):
                    problem.setdefault('label_shape', []).append(lab)
        if problem:
            bad.append({'identity': _result_identity(config, r), 'problem': problem})
    checks.add('labels.results', not bad, 'every returned context with its labels', bad[:10])
    by_label = entry.get('by_label')
    want_by = [{'column': b['column'], 'contexts': b['contexts'],
                'contexts_suffix': b['contexts_suffix'], 'occurrences': b['occurrences']}
               for b in view['by_label']]
    complete = released is not None and len(released) == len(exp.ctxs) and not stopped
    got_by = None
    if isinstance(by_label, list):
        got_by = []
        for b in by_label:
            occ = b.get('occurrences') or {}
            got_by.append({'column': b.get('column'),
                           'contexts': (b.get('contexts') or {}).get('value'),
                           'contexts_suffix': (b.get('contexts_suffix') or {}).get('value'),
                           'occurrences': occ.get('value')})
    checks.add('labels.by_label', got_by == want_by, want_by, got_by)
    if complete:
        judge_count(checks, 'labels.counts.labels', counts.get('labels'), view['labels'], True,
                    unit='labels')
        if config.placement == 'record':
            judge_count(checks, 'labels.counts.occurrences', counts.get('occurrences'),
                        view['occurrences'], True, unit='placed_occurrences')
    elif released is not None and not stopped:
        # a list cut at max_contexts: the labels of the returned contexts, at_least (SPEC 14.6)
        lab = counts.get('labels') or {}
        checks.add('labels.counts.labels', lab.get('relation') == 'at_least'
                   and lab.get('value') == view['labels'],
                   {'relation': 'at_least', 'value': view['labels']}, lab)
        if config.placement == 'record':
            occ = counts.get('occurrences') or {}
            checks.add('labels.counts.occurrences', occ.get('relation') == 'at_least'
                       and occ.get('value') == view['occurrences'],
                       {'relation': 'at_least', 'value': view['occurrences']}, occ)
    else:
        checks.mark('labels.counts', INCONCLUSIVE, note='a stopped list: label counts not '
                    'judged as exact')


# -- the whole /pattern answer --------------------------------------------------------------

def sticky_stop(stop, request):
    """The stop every later pattern of the request carries after |stop| (SPEC section 7.6: a
    budget stop is sticky, phase discovery; section 21.6: a request-wide stop before the pattern
    names level 0); None for a stop that ends only its own pattern (stop_at_threshold)."""
    if not isinstance(stop, dict) or stop.get('reason') not in BUDGET_REASONS:
        return None
    out = {'phase': 'discovery', 'reason': stop['reason']}
    if request.get('max_mismatches', 0):
        out['mismatch_level'] = 0
    return out


def compare_pattern(index, config, request, answer, max_contexts=None):
    """-> [(position, Checks, Expectation)] for every pattern of the request."""
    out = []
    patterns = answer.get('patterns') if isinstance(answer, dict) else None
    if not isinstance(patterns, list) or len(patterns) != len(request['patterns']):
        c = Checks()
        c.add('answer.patterns', False, len(request['patterns']),
              None if not isinstance(patterns, list) else len(patterns))
        return [(None, c, None)]
    sticky = None
    steps = 0
    for i, entry in enumerate(patterns):
        checks, exp = compare_entry(index, config, request, i, entry, sticky, max_contexts)
        stop = entry.get('stop') if isinstance(entry, dict) else None
        if sticky is None:
            sticky = sticky_stop(stop, request)
        work = entry.get('work') if isinstance(entry, dict) else None
        if isinstance(work, dict) and isinstance(work.get('steps'), int):
            steps += work['steps']
            sbm = work.get('steps_by_mismatches')
            if isinstance(sbm, dict):
                checks.add('stop.steps_by_mismatches_sum', sum(sbm.values()) == work['steps'],
                           work['steps'], sbm)
        out.append((i, checks, exp))
    if 'max_steps' in request and out:
        out[-1][1].add('stop.steps_within_budget', steps <= request['max_steps'],
                       f'<= {request["max_steps"]}', steps)
    return out


# -- /resolve -----------------------------------------------------------------------------

def compare_resolve(index, config, request, answer, status=200):
    """Checks of a /resolve answer (|status| its HTTP status)."""
    checks = Checks()
    support = request.get('support', 'kmer')
    try:
        want = O.resolve(index, request['sequence'], request['labels'],
                         request.get('min_block_kmers', 1), support)
    except ValueError as e:
        # a name that is no column and no indexed record (a record shorter than k is not
        # indexed): the request is refused, naming the label
        name = str(e).split("'")[1] if "'" in str(e) else ''
        text = answer.get('error', '') if isinstance(answer, dict) else ''
        checks.add('resolve.unknown_label', status == 400 and name in text,
                   f'400 naming {name!r}', [status, answer])
        return checks, None
    if status != 200:
        checks.add('resolve.status', False, 200, [status, answer])
        return checks, want
    checks.add('resolve.k', answer.get('k') == want['k'], want['k'], answer.get('k'))
    checks.add('resolve.regime', answer.get('regime') == config.mode, config.mode,
               answer.get('regime'))
    checks.add('resolve.support', answer.get('support') == support, support, answer.get('support'))
    checks.add('resolve.num_kmers', answer.get('num_kmers') == want['num_kmers'],
               want['num_kmers'], answer.get('num_kmers'))
    checks.add('resolve.graph_runs', answer.get('graph_runs') == want['graph_runs'],
               want['graph_runs'], answer.get('graph_runs'))
    got_labels = answer.get('labels')
    if isinstance(got_labels, list):
        got = [{f: l.get(f) for f in ('label', 'kind', 'kmers_supported', 'runs')}
               for l in got_labels]
    else:
        got = got_labels
    checks.add('resolve.labels', got == want['labels'], want['labels'], got)
    cands = answer.get('candidates')
    if isinstance(cands, list):
        gotc = [{f: c.get(f) for f in ('kmer_interval', 'bp_interval', 'labels')} for c in cands]
    else:
        gotc = cands
    checks.add('resolve.candidates', gotc == want['candidates'], want['candidates'], gotc)
    checks.add('resolve.labels_truncated', answer.get('labels_truncated') is None, None,
               answer.get('labels_truncated'))
    return checks, want
