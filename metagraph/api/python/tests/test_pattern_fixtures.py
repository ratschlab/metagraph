"""The POST /pattern fixture bodies (data/traverse/pattern/, written by
scripts/traversal/pattern_fixtures.py) against docs/SPEC-pattern-search.md, contract version 1
(milestone 1, and increment 3's output.labels "all", SPEC §14):

  - the field lists: every table of the SPEC marked `<!-- schema: NAME -->` names exactly the
    fields SCHEMA[NAME] below knows, so that the SPEC and this check cannot drift apart;
  - every answer: exactly the fields its shape has (answered entry, error slot, refusal,
    capabilities block), the types and enumerations of the SPEC, and the rules that tie the
    fields together (a count's relation and its value, the offsets of the scope, the order of
    the results, what withheld, cut and retrieval_complete imply) -- checked on the bodies,
    never by asking a server;
  - the fixture set: it covers every mode, scope, strand setting, withheld and cut reason, slot
    error and refusal code a service has to handle, and the README names each fixture.

No metagraph import and no server: the service copies these bodies into its own tests, and
this check is what makes them a contract rather than a sample.
"""

import json
import math
import os
import re
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data', 'traverse', 'pattern')
REPO = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
SPEC = os.path.join(REPO, 'docs', 'SPEC-pattern-search.md')
PATTERN_CPP = os.path.join(REPO, 'src', 'cli', 'pattern.cpp')
SERVER_UTILS_CPP = os.path.join(REPO, 'src', 'cli', 'server_utils.cpp')
PROMPT = os.path.join(REPO, 'docs', 'PROMPT-search-service-pattern.md')

IUPAC = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'R': 'AG', 'Y': 'CT', 'S': 'CG', 'W': 'AT',
         'K': 'GT', 'M': 'AC', 'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'}
COMPLEMENT = str.maketrans('ACGTRYSWKMBDHVN', 'TGCAYRSWMKVHDBN')

# the enumerations of the SPEC (§§6-10)
MODES = ('count', 'all_or_count', 'partial')
SCOPES = ('suffix', 'any_offset', 'long')
ABSENCE_SCOPES = {'suffix': 'suffix_only', 'any_offset': 'any_offset', 'long': 'long'}
RELATIONS = ('exact', 'at_least', 'bounds', 'unknown')
UNITS = ('graph_contexts', 'anchors', 'paths', 'placed_occurrences', 'labels')
GRAPH_MODES = ('basic', 'canonical', 'primary')
STRAND_SYMBOLS = {'forward': '+', 'reverse': '-', 'palindromic': '='}
STRAND_KEYS = {'+': '+', '-': '-', '=': 'both'}
# increment 3 (output.labels "all", SPEC §14) adds the annotation phases and reasons
ANNOTATION_PHASES = ('label_discovery', 'placement', 'output')
STOP_PHASES = ('discovery', 'mask_scan', 'extraction') + ANNOTATION_PHASES
STOP_REASONS = ('max_steps', 'time', 'max_contexts', 'max_anchors', 'max_annotation_work',
                'max_memory')
WITHHELD = ('count_above_threshold', 'threshold_crossed', 'discovery_budget', 'deadline',
            'paths_later_increment', 'annotation_budget', 'anchor_labels_truncated',
            'output_budget')
CUT = ('max_contexts', 'max_steps', 'time', 'max_anchors', 'max_memory')
NOTES = ('low_complexity_pattern', 'strand_unknown_canonical', 'paths_later_increment',
         'annotation_unbudgeted', 'record_bounds_unknown', 'annotation_not_read')
SLOT_ERRORS = ('bad_alphabet', 'information_below_floor', 'scope_unsupported')
REFUSALS = ('invalid_request', 'later_increment', 'resident_only', 'mask_required',
            'representation_unsupported', 'primary_unwrapped', 'alphabet_unsupported', 'deadline',
            'annotation_unbudgeted')
UNAVAILABLE = ('mask_required', 'representation_unsupported', 'primary_unwrapped',
               'alphabet_unsupported', 'multi_graph_later_increment')
MASKS = ('file', 'built_at_load', 'absent')
PLACEMENTS = ('record', 'global', 'none', 'none_canonical')
# an entry's placement (increment 3): the index's, or not_requested (output.occurrences false)
ENTRY_PLACEMENTS = PLACEMENTS + ('not_requested',)
LABELS_STATUS = ('complete', 'truncated', 'refused', 'not_read', 'output_budget')
SUPPORTS = ('record_verified', 'label_intersection')
ANNOTATIONS = ('budgeted', 'unbudgeted')
ORIENTATION_RANK = {'+': 0, '-': 1, '=': 2, 'forward': 0, 'reverse': 1, 'palindromic': 2}

# The fields of each object of the SPEC, by the name of its schema table. The checks below
# decide when each is present; this is the list the SPEC's tables must name exactly.
SCHEMA = {
    'request': ['patterns', 'mode', 'scope', 'strands', 'stop_at_threshold', 'max_contexts',
                'max_anchors', 'max_steps', 'time_budget_ms', 'output',
                'max_labels_per_anchor', 'max_annotation_work', 'max_memory_mb', 'max_labels',
                'max_occurrences_per_label', 'allow_unbudgeted_annotation'],
    'request_pattern': ['id', 'dna', 'iupac'],
    'request_output': ['labels', 'occurrences', 'paths'],
    'refusal': ['error', 'code'],
    'answer': ['pattern_contract_version', 'mode', 'output', 'index', 'limits', 'timing',
               'patterns'],
    'index': ['index_ns', 'index_fp', 'release', 'k', 'graph_mode', 'alphabet',
              'strand_stated'],
    'limits': ['max_contexts', 'max_anchors', 'max_steps', 'time_budget_ms',
               'finalize_reserve_ms', 'min_information_bits', 'max_patterns',
               'stop_at_threshold', 'max_labels_per_anchor', 'max_annotation_work',
               'max_memory_mb', 'max_labels', 'max_occurrences_per_label',
               'allow_unbudgeted_annotation', 'clamped'],
    'clamped': ['field', 'requested', 'effective'],
    'timing': ['elapsed_ms', 'label_discovery_ms', 'placement_ms'],
    'entry': ['id', 'kind', 'pattern', 'length', 'information_bits', 'anchor_information_bits',
              'error', 'mode', 'scope', 'strands', 'palindromic', 'counts', 'work', 'stop',
              'retrieval_complete', 'withheld', 'returned', 'cut', 'results', 'absence_scope',
              'determinism', 'notes', 'timing', 'placement', 'annotation', 'by_label',
              'rows_refused', 'anchors_truncated', 'labels_cut', 'occurrences_cut'],
    'counts': ['contexts', 'anchors', 'paths', 'labels', 'occurrences'],
    'count': ['value', 'relation', 'unit', 'lower', 'upper'],
    'contexts_count': ['suffix', 'by_offset', 'by_strand', 'by_orientation'],
    'anchors_count': ['by_strand', 'by_orientation'],
    'work': ['ranges_visited', 'mask_scans', 'steps', 'annotation_rows', 'annotation_units',
             'memory_bytes'],
    'stop': ['phase', 'reason'],
    'reason': ['reason'],
    'result': ['kmer', 'instance', 'offset', 'strand', 'orientation', 'node', 'row', 'support',
               'labels_status', 'labels_total', 'labels'],
    'error': ['code', 'message'],
    'capabilities': ['pattern_contract_version', 'available', 'unavailable_reason', 'modes',
                     'default_mode', 'projections', 'default_projection',
                     'projections_later_increment', 'kinds', 'kinds_later_increment',
                     'default_scope', 'scopes_by_graph_mode', 'scopes', 'long_patterns',
                     'strands', 'default_strands', 'graph_cleaned', 'records_shorter_than_k',
                     'resident_only', 'caps', 'default_time_budget_ms', 'finalize_reserve_ms',
                     'caps_rule', 'graph_mode', 'k', 'alphabet', 'strand_stated', 'mask',
                     'placement', 'support', 'annotation', 'default_occurrences'],
    'capabilities_multi': ['pattern_contract_version', 'available', 'unavailable_reason'],
    # increment 3 (SPEC §14)
    'anchor_truncated': ['kmer', 'row', 'cap', 'total'],
    'occurrence': ['seq_id', 'record', 'strand', 'nt_coords', 'nt_length'],
    'occurrence_global': ['kmer_coord', 'offset', 'strand'],
    'row_refused': ['kmer', 'row', 'phase', 'reason', 'needed_bytes', 'available_bytes'],
    'label': ['column', 'support', 'occurrences', 'occurrence_list'],
    'by_label': ['graph', 'column', 'contexts', 'contexts_suffix', 'occurrences'],
    'labels_cut': ['reason', 'returned'],
    'occurrences_cut': ['reason', 'labels'],
}
ENTRY_DESCRIPTION = ['id', 'kind', 'pattern', 'length', 'information_bits',
                     'anchor_information_bits']
ENTRY_ANSWERED = ENTRY_DESCRIPTION + ['mode', 'scope', 'strands', 'palindromic', 'counts',
                                      'work', 'stop', 'retrieval_complete', 'absence_scope',
                                      'determinism', 'notes', 'timing']
ENTRY_RETRIEVAL = ['withheld', 'returned', 'cut', 'results']
# output.labels "all" in a retrieval mode (increment 3)
ENTRY_LABELS = ['placement', 'annotation', 'by_label', 'rows_refused', 'anchors_truncated',
                'labels_cut', 'occurrences_cut']
LIMITS_LABELS = ['max_labels_per_anchor', 'max_annotation_work', 'max_memory_mb', 'max_labels',
                 'max_occurrences_per_label']
CAPS = ['max_contexts', 'max_anchors', 'max_steps', 'time_budget_ms', 'min_information_bits',
        'max_patterns'] + LIMITS_LABELS


def load(*path):
    with open(os.path.join(DATA, *path), encoding='utf-8') as f:
        return json.load(f)


def spec_tables(text):
    """{name: [field, ...]}: the first column of the table after each schema marker."""
    out = {}
    lines = text.split('\n')
    for i, line in enumerate(lines):
        m = re.match(r'<!-- schema: (\w+) -->\s*$', line.strip())
        if not m:
            continue
        rows = []
        j = i + 1
        while j < len(lines) and not lines[j].startswith('|'):
            j += 1
        for row in lines[j + 2:]:       # past the header and the separator
            if not row.startswith('|'):
                break
            first = row.split('|')[1].strip()
            name = re.fullmatch(r'`([^`]+)`', first)
            assert name, f'schema {m.group(1)}: first column {first!r} is not a `field`'
            rows.append(name.group(1))
        assert m.group(1) not in out, f'schema {m.group(1)} twice'
        out[m.group(1)] = rows
    return out


def is_int(x):
    return isinstance(x, int) and not isinstance(x, bool)


def is_num(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def information_bits(pattern):
    return sum(math.log2(4 / len(IUPAC[c])) for c in pattern)


def matches(pattern, text):
    return len(pattern) == len(text) and all(t in IUPAC[c] for c, t in zip(pattern, text))


def revcomp(s):
    return s.translate(COMPLEMENT)[::-1]


class Checker:
    """Validates one answer body; every failure names the fixture and the JSON path."""

    def __init__(self, test, name):
        self.t = test
        self.name = name

    def ok(self, condition, path, what=''):
        if not condition:
            self.t.fail(f'{self.name}: {path}: {what}')

    def keys(self, obj, expected, path):
        self.ok(isinstance(obj, dict), path, f'expected an object, got {obj!r}')
        self.ok(set(obj) == set(expected), path,
                f'fields {sorted(obj)}, expected {sorted(expected)}')

    def one_of(self, value, values, path):
        self.ok(value in values, path, f'{value!r} not in {values}')

    # ------------------------------------------------------------- counts

    def count(self, c, unit, path, extra=()):
        keys = ['value', 'relation', 'unit'] + list(extra)
        if isinstance(c, dict) and c.get('relation') == 'bounds':
            keys += ['lower', 'upper']
        self.keys(c, keys, path)
        self.one_of(c['relation'], RELATIONS, path + '.relation')
        self.ok(c['unit'] == unit, path + '.unit', f'{c["unit"]!r}, expected {unit!r}')
        if c['relation'] == 'unknown':
            self.ok(c['value'] is None, path + '.value', 'unknown counts are null')
        else:
            self.ok(is_int(c['value']) and c['value'] >= 0, path + '.value', repr(c['value']))
        if c['relation'] == 'bounds':
            self.ok(is_int(c['lower']) and is_int(c['upper']) and c['lower'] <= c['upper']
                    and c['value'] == c['lower'], path, 'bounds: value == lower <= upper')

    def summed(self, total, parts, path):
        """A total over disjoint parts: exact iff every part is, and then their sum."""
        if total['relation'] == 'exact':
            self.ok(all(p['relation'] == 'exact' for p in parts), path,
                    'an exact total has exact parts')
            self.ok(total['value'] == sum(p['value'] for p in parts), path,
                    'an exact total is the sum of its parts')
        if all(p['relation'] == 'exact' for p in parts):
            self.ok(total['relation'] == 'exact', path, 'exact parts make an exact total')
        if total['relation'] == 'unknown':
            self.ok(all(p['relation'] == 'unknown' for p in parts), path,
                    'an unknown total has unknown parts')

    def by_orientation(self, c, entry, strand_stated, unit, path):
        key = 'by_strand' if strand_stated else 'by_orientation'
        other = 'by_orientation' if strand_stated else 'by_strand'
        self.ok(other not in c, path, f'{other} on a graph with strand_stated={strand_stated}')
        expected = [STRAND_KEYS[s] for s in entry['strands']] if strand_stated \
            else list(entry['strands'])
        self.keys(c[key], expected, f'{path}.{key}')
        for k, v in c[key].items():
            self.count(v, unit, f'{path}.{key}.{k}')
        self.summed(c, list(c[key].values()), f'{path}.{key}')

    # ------------------------------------------------------------- answers

    def answer(self, a, request, capabilities=None):
        self.keys(a, SCHEMA['answer'], 'answer')
        self.ok(a['pattern_contract_version'] == 1, 'pattern_contract_version')
        mode = request.get('mode', 'all_or_count')
        self.ok(a['mode'] == mode, 'mode', f'{a["mode"]!r}, the request asks {mode!r}')
        asked = request.get('output', {})
        labelled = mode != 'count' and asked.get('labels') == 'all'
        if mode == 'count':
            self.ok(a['output'] is None, 'output', 'null in mode count')
        elif labelled:
            self.ok(a['output'] == {'labels': 'all',
                                    'occurrences': asked.get('occurrences', True)},
                    'output', repr(a['output']))
        else:
            self.ok(a['output'] == {'labels': 'none'}, 'output', repr(a['output']))

        index = a['index']
        self.keys(index, SCHEMA['index'], 'index')
        for f in ('index_ns', 'index_fp'):
            self.ok(index[f] is None or isinstance(index[f], str), 'index.' + f)
        self.ok(isinstance(index['release'], str), 'index.release')
        self.ok(is_int(index['k']) and index['k'] >= 2, 'index.k')
        self.one_of(index['graph_mode'], GRAPH_MODES, 'index.graph_mode')
        self.one_of(index['alphabet'], ('$ACGT', '$ACGTN'), 'index.alphabet')
        self.ok(index['strand_stated'] is (index['graph_mode'] == 'basic'),
                'index.strand_stated', 'true exactly on basic graphs')

        limits = a['limits']
        # the annotation limits only in the answers that read annotation (increment 3)
        self.keys(limits, [f for f in SCHEMA['limits']
                           if labelled or f not in LIMITS_LABELS + ['allow_unbudgeted_annotation']],
                  'limits')
        if labelled:
            for f in LIMITS_LABELS:
                self.ok(is_int(limits[f]) and limits[f] >= 0, 'limits.' + f)
            self.ok(limits['allow_unbudgeted_annotation']
                    is request.get('allow_unbudgeted_annotation', False),
                    'limits.allow_unbudgeted_annotation')
        for f in ('max_contexts', 'max_anchors', 'max_steps', 'max_patterns'):
            self.ok(is_int(limits[f]) and limits[f] >= 0, 'limits.' + f)
        self.ok(limits['max_steps'] >= 1, 'limits.max_steps')
        for f in ('time_budget_ms', 'finalize_reserve_ms', 'min_information_bits'):
            self.ok(is_num(limits[f]), 'limits.' + f)
        self.ok(limits['time_budget_ms'] > limits['finalize_reserve_ms'], 'limits.time_budget_ms',
                'above the finalisation reserve')
        self.ok(limits['stop_at_threshold'] is request.get('stop_at_threshold', False),
                'limits.stop_at_threshold')
        self.ok(isinstance(limits['clamped'], list), 'limits.clamped')
        order = ['max_contexts', 'max_anchors', 'max_steps', 'time_budget_ms'] + LIMITS_LABELS
        clamped = []
        for i, c in enumerate(limits['clamped']):
            path = f'limits.clamped[{i}]'
            self.keys(c, SCHEMA['clamped'], path)
            self.one_of(c['field'], order, path + '.field')
            self.ok(is_num(c['requested']) and is_num(c['effective'])
                    and c['requested'] > c['effective'], path, 'requested above effective')
            self.ok(request.get(c['field']) == c['requested'], path, 'the value asked for')
            if c['field'] in limits:
                self.ok(limits[c['field']] == c['effective'], path, 'the effective limit')
            clamped.append(c['field'])
        self.ok(clamped == [f for f in order if f in clamped], 'limits.clamped', 'order')
        echoed = [f for f in order if f in limits]
        for f in echoed:
            if f in request and f not in clamped:
                self.ok(limits[f] == request[f], 'limits.' + f, 'a value under its cap is kept')
        if capabilities is not None:
            caps = capabilities['caps']
            for f in echoed:
                if f not in request:
                    expected = capabilities['default_time_budget_ms'] \
                        if f == 'time_budget_ms' else caps[f]
                    self.ok(limits[f] == expected, 'limits.' + f, 'the default')
                self.ok(limits[f] <= caps[f], 'limits.' + f, 'within the cap')
            self.ok(limits['finalize_reserve_ms'] == capabilities['finalize_reserve_ms'],
                    'limits.finalize_reserve_ms')
            self.ok(limits['min_information_bits'] == caps['min_information_bits'],
                    'limits.min_information_bits')
            self.ok(limits['max_patterns'] == caps['max_patterns'], 'limits.max_patterns')

        self.timing(a['timing'], 'timing')
        self.ok(len(a['patterns']) == len(request['patterns']), 'patterns',
                'one entry per request pattern')
        for i, (entry, asked) in enumerate(zip(a['patterns'], request['patterns'])):
            self.entry(entry, asked, request, a, f'patterns[{i}]')

    def timing(self, t, path, labelled=False):
        self.keys(t, SCHEMA['timing'] if labelled else ['elapsed_ms'], path)
        for k, v in t.items():
            self.ok(is_num(v) and v >= 0, f'{path}.{k}')

    def entry(self, e, asked, request, a, path):
        k = a['index']['k']
        strand_stated = a['index']['strand_stated']
        mode = a['mode']
        labelled = isinstance(a['output'], dict) and a['output']['labels'] == 'all'
        # the request named the labels (or an annotation field) and this answer reads none
        named = (request.get('output', {}).get('labels') == 'all'
                 or any(f in request for f in LIMITS_LABELS + ['allow_unbudgeted_annotation']))
        self.ok(set(e) <= set(SCHEMA['entry']), path, f'fields outside the SPEC: '
                f'{sorted(set(e) - set(SCHEMA["entry"]))}')
        self.ok(e['id'] == asked.get('id'), path + '.id', 'echoed')
        kind = 'dna' if 'dna' in asked else 'iupac'
        self.ok(e['kind'] == kind, path + '.kind')
        text = asked[kind].upper()

        if 'error' in e:
            err = e['error']
            self.keys(err, SCHEMA['error'], path + '.error')
            self.one_of(err['code'], SLOT_ERRORS, path + '.error.code')
            self.ok(isinstance(err['message'], str) and err['message'], path + '.error.message')
            if err['code'] == 'bad_alphabet':
                self.keys(e, ['id', 'kind', 'error'], path)
                alphabet = 'ACGT' if kind == 'dna' else ''.join(IUPAC)
                self.ok(not text or any(c not in alphabet for c in text), path,
                        'bad_alphabet for a pattern inside the alphabet')
                return
            self.keys(e, ENTRY_DESCRIPTION + ['error'], path)
            self.description(e, text, k, path)
            if err['code'] == 'information_below_floor':
                bits = e['anchor_information_bits'] if e['length'] > k else e['information_bits']
                self.ok(bits < a['limits']['min_information_bits'], path,
                        'refused above the floor')
            return

        expected = ENTRY_ANSWERED + (ENTRY_RETRIEVAL if mode != 'count' else []) \
            + (ENTRY_LABELS if labelled else [])
        self.keys(e, expected, path)
        self.description(e, text, k, path)
        self.ok(e['mode'] == mode, path + '.mode')
        L = e['length']
        scope = 'long' if L > k else request.get('scope', 'any_offset')
        self.ok(e['scope'] == scope, path + '.scope', f'{e["scope"]!r}, expected {scope!r}')
        self.ok(e['absence_scope'] == ABSENCE_SCOPES[scope], path + '.absence_scope')
        palindromic = revcomp(text) == text
        self.ok(e['palindromic'] is palindromic, path + '.palindromic')
        strands = request.get('strands', 'both')
        if palindromic:
            orientations = ['palindromic']
        else:
            orientations = {'both': ['forward', 'reverse'], 'forward': ['forward'],
                            'reverse': ['reverse']}[strands]
        expected_strands = [STRAND_SYMBOLS[o] for o in orientations] if strand_stated \
            else orientations
        self.ok(e['strands'] == expected_strands, path + '.strands',
                f'{e["strands"]}, expected {expected_strands}')

        # counts
        counts = e['counts']
        if L > k:
            self.keys(counts, ['anchors', 'paths', 'labels', 'occurrences'], path + '.counts')
            anchors = counts['anchors']
            self.count(anchors, 'anchors', path + '.counts.anchors',
                       extra=['by_strand' if strand_stated else 'by_orientation'])
            self.by_orientation(anchors, e, strand_stated, 'anchors', path + '.counts.anchors')
            self.count(counts['paths'], 'paths', path + '.counts.paths')
            if (anchors['relation'], anchors['value']) == ('exact', 0):
                self.ok(counts['paths'] == {'value': 0, 'relation': 'exact', 'unit': 'paths'},
                        path + '.counts.paths', 'no anchor, no path: exact 0')
            else:
                self.ok(counts['paths']['relation'] == 'unknown', path + '.counts.paths',
                        'paths are a later increment')
            total = anchors
        else:
            self.keys(counts, ['contexts', 'labels', 'occurrences'], path + '.counts')
            c = counts['contexts']
            self.count(c, 'graph_contexts', path + '.counts.contexts',
                       extra=['suffix', 'by_offset',
                              'by_strand' if strand_stated else 'by_orientation'])
            self.count(c['suffix'], 'graph_contexts', path + '.counts.contexts.suffix')
            offsets = [k - L] if scope == 'suffix' else list(range(k - L + 1))
            self.keys(c['by_offset'], [str(o) for o in offsets],
                      path + '.counts.contexts.by_offset')
            for o, v in c['by_offset'].items():
                self.count(v, 'graph_contexts', f'{path}.counts.contexts.by_offset.{o}')
            self.ok(c['suffix'] == c['by_offset'][str(k - L)], path + '.counts.contexts.suffix',
                    'the count at offset k - L')
            self.summed(c, list(c['by_offset'].values()), path + '.counts.contexts.by_offset')
            self.by_orientation(c, e, strand_stated, 'graph_contexts', path + '.counts.contexts')
            total = c
        for f, unit in (('labels', 'labels'), ('occurrences', 'placed_occurrences')):
            if labelled:
                self.count(counts[f], unit, f'{path}.counts.{f}')
            else:
                self.ok(counts[f] == {'value': None, 'relation': 'unknown', 'unit': unit},
                        f'{path}.counts.{f}', 'no annotation is read without labels "all"')

        # work, stop, determinism, notes
        w = e['work']
        self.keys(w, SCHEMA['work'] if labelled else ['ranges_visited', 'mask_scans', 'steps'],
                  path + '.work')
        self.ok(all(is_int(v) and v >= 0 for v in w.values()), path + '.work')
        self.ok(w['steps'] >= w['ranges_visited'], path + '.work', 'steps >= ranges_visited')
        stop = e['stop']
        if stop is not None:
            self.keys(stop, SCHEMA['stop'], path + '.stop')
            self.one_of(stop['phase'], STOP_PHASES, path + '.stop.phase')
            self.one_of(stop['reason'], STOP_REASONS, path + '.stop.reason')
            self.ok(total['relation'] != 'exact'
                    or stop['phase'] in ('extraction',) + ANNOTATION_PHASES,
                    path + '.stop', 'a stop in discovery leaves no exact count')
            self.ok(stop['phase'] not in ANNOTATION_PHASES or labelled, path + '.stop.phase',
                    'an annotation phase without labels "all"')
            if stop['reason'] == 'max_steps':
                self.ok(w['steps'] <= a['limits']['max_steps'], path + '.work.steps')
        else:
            self.ok(total['relation'] == 'exact', path + '.counts',
                    'a count without a stop is exact')
        self.one_of(e['determinism'], ('full', 'time_limited'), path + '.determinism')
        if stop is not None and stop['reason'] == 'time':
            self.ok(e['determinism'] == 'time_limited', path + '.determinism',
                    'a time stop is not deterministic')
        if e['determinism'] == 'time_limited':
            # (with labels "all", a time stop of the reads after a step stop of the engine:
            # the engine's stop is the one stated)
            self.ok(stop is not None and (stop['reason'] == 'time' or labelled), path + '.stop')
        self.ok(isinstance(e['notes'], list) and all(n in NOTES for n in e['notes']),
                path + '.notes', repr(e['notes']))
        self.ok(e['notes'] == [n for n in NOTES if n in e['notes']], path + '.notes', 'order')
        self.ok(('strand_unknown_canonical' in e['notes']) is (not strand_stated),
                path + '.notes', 'strand_unknown_canonical exactly where no strand is known')
        self.ok(('paths_later_increment' in e['notes']) is (L > k), path + '.notes')
        self.ok(('annotation_not_read' in e['notes']) is (named and not labelled),
                path + '.notes', 'annotation_not_read exactly where labels were named, not read')
        if labelled:
            self.ok(('annotation_unbudgeted' in e['notes'])
                    is (e['annotation'] == 'unbudgeted'), path + '.notes')
            self.ok(('record_bounds_unknown' in e['notes']) is (e['placement'] == 'global'),
                    path + '.notes')
        else:
            self.ok(not {'annotation_unbudgeted', 'record_bounds_unknown'} & set(e['notes']),
                    path + '.notes')
        self.timing(e['timing'], path + '.timing', labelled)

        if mode == 'count':
            self.ok(e['retrieval_complete'] is False, path + '.retrieval_complete',
                    'a count returns no context')
            return
        self.retrieval(e, text, total, request, a, path)

    def description(self, e, text, k, path):
        self.ok(e['pattern'] == text, path + '.pattern', 'the pattern in upper case')
        self.ok(e['length'] == len(text), path + '.length')
        self.ok(is_num(e['information_bits'])
                and abs(e['information_bits'] - information_bits(text)) < 1e-9,
                path + '.information_bits')
        if len(text) > k:
            self.ok(is_num(e['anchor_information_bits'])
                    and abs(e['anchor_information_bits'] - information_bits(text[:k])) < 1e-9,
                    path + '.anchor_information_bits', 'the bits of [0, k)')
        else:
            self.ok(e['anchor_information_bits'] is None, path + '.anchor_information_bits')

    def retrieval(self, e, text, total, request, a, path):
        k = a['index']['k']
        strand_stated = a['index']['strand_stated']
        mode = a['mode']
        labelled = isinstance(a['output'], dict) and a['output']['labels'] == 'all'
        L = len(text)
        results = e['results']
        self.ok(isinstance(results, list) and e['returned'] == len(results), path + '.returned',
                'the length of results')
        for f in ('withheld', 'cut'):
            if e[f] is not None:
                self.keys(e[f], SCHEMA['reason'], f'{path}.{f}')
        if e['withheld'] is not None:
            self.one_of(e['withheld']['reason'], WITHHELD, path + '.withheld.reason')
            self.ok(results == [] and e['cut'] is None and not e['retrieval_complete'],
                    path, 'withheld: no results, no cut, not complete')
        if e['cut'] is not None:
            self.one_of(e['cut']['reason'], CUT, path + '.cut.reason')
            self.ok(mode == 'partial' and not e['retrieval_complete'], path + '.cut',
                    'only partial cuts, and a cut list is not complete')
        if e['retrieval_complete']:
            self.ok(total['relation'] == 'exact' and e['stop'] is None
                    and e['withheld'] is None and e['cut'] is None, path,
                    'complete: exact, no stop, nothing withheld or cut')
            self.ok(e['returned'] == (total['value'] if L <= k else 0), path,
                    'complete: every context returned')
            if L > k:
                self.ok(total['value'] == 0, path, 'a long pattern is complete only without '
                        'anchors')
        elif e['withheld'] is None and mode == 'partial':
            labels_stated = labelled and (
                e['rows_refused'] or e['anchors_truncated'] or e['stop'] is not None
                or e['labels_cut'] is not None or e['occurrences_cut'] is not None)
            self.ok(e['cut'] is not None or labels_stated, path + '.cut',
                    'an incomplete partial list states its cut')
        if mode == 'all_or_count':
            self.ok(e['retrieval_complete'] or e['withheld'] is not None, path,
                    'all_or_count returns all, or withholds')
        if L > k:
            self.ok(results == [], path + '.results', 'no result for a long pattern in version 1')
            if not e['retrieval_complete']:
                self.ok(e['withheld'] == {'reason': 'paths_later_increment'}, path + '.withheld')
        if mode == 'partial':
            self.ok(e['returned'] <= a['limits']['max_contexts'], path + '.returned', 'the cap')

        key = 'strand' if strand_stated else 'orientation'
        allowed = set(e['strands'])
        order = []
        for i, r in enumerate(results):
            rp = f'{path}.results[{i}]'
            self.keys(r, ['kmer', 'instance', 'offset', key, 'node', 'row']
                      + (['support', 'labels_status', 'labels_total', 'labels'] if labelled
                         else []), rp)
            self.ok(isinstance(r['kmer'], str) and len(r['kmer']) == k
                    and set(r['kmer']) <= set('ACGTN'), rp + '.kmer')
            self.ok(is_int(r['offset']) and 0 <= r['offset'] <= k - L, rp + '.offset')
            if e['scope'] == 'suffix':
                self.ok(r['offset'] == k - L, rp + '.offset', 'suffix scope: offset k - L')
            self.ok(r['instance'] == r['kmer'][r['offset']:r['offset'] + L], rp + '.instance')
            self.ok(r[key] in allowed, rp + '.' + key, f'{r[key]!r} not searched')
            oriented = revcomp(text) if r[key] in ('-', 'reverse') else text
            self.ok(matches(oriented, r['instance']), rp + '.instance',
                    f'{r["instance"]} does not instantiate {oriented}')
            self.ok(is_int(r['node']) and r['node'] >= 1, rp + '.node')
            self.ok(r['row'] is None or (is_int(r['row']) and r['row'] >= 0), rp + '.row')
            if a['index']['graph_mode'] == 'basic' and r['row'] is not None:
                self.ok(r['row'] == r['node'] - 1, rp + '.row', 'node - 1 on a basic graph')
            order.append((r['node'], r['offset'], ORIENTATION_RANK[r[key]]))
        self.ok(order == sorted(order) and len(set(order)) == len(order), path + '.results',
                'ordered by (node, offset, orientation), each context once')
        if labelled:
            self.labels(e, total, L, a, path)

    def labels(self, e, total, L, a, path):
        """output.labels "all" (SPEC §14): the label fields, their values and the rules that tie
        them to retrieval_complete, withheld and the counts."""
        k = a['index']['k']
        mode = a['mode']
        self.one_of(e['placement'], ENTRY_PLACEMENTS, path + '.placement')
        if a['index']['graph_mode'] != 'basic':
            self.ok(e['placement'] in ('none_canonical', 'not_requested'), path + '.placement')
        self.one_of(e['annotation'], ANNOTATIONS, path + '.annotation')
        records = e['placement'] == 'record'
        placed = e['placement'] in ('record', 'global')
        for i, t in enumerate(e['anchors_truncated']):
            tp = f'{path}.anchors_truncated[{i}]'
            self.keys(t, SCHEMA['anchor_truncated'], tp)
            self.ok(t['cap'] == a['limits']['max_labels_per_anchor'] and t['total'] > t['cap'],
                    tp, 'cut at the cap')
        for i, x in enumerate(e['rows_refused']):
            xp = f'{path}.rows_refused[{i}]'
            self.keys(x, SCHEMA['row_refused'], xp)
            self.one_of(x['phase'], ('label_discovery', 'placement'), xp + '.phase')
            self.ok(x['reason'] == 'max_memory', xp + '.reason')
        for f, schema in (('labels_cut', 'labels_cut'), ('occurrences_cut', 'occurrences_cut')):
            if e[f] is not None:
                self.keys(e[f], SCHEMA[schema], f'{path}.{f}')
                self.ok(mode == 'partial' and not e['retrieval_complete'], f'{path}.{f}',
                        'only partial cuts its lists, and a cut answer is not complete')
        if e['withheld'] is not None:
            self.ok(e['by_label'] is None, path + '.by_label', 'null when withheld')
            for f in ('labels', 'occurrences'):
                self.ok(e['counts'][f]['relation'] == 'unknown', f'{path}.counts.{f}',
                        'nothing published when withheld')
            return
        self.ok(isinstance(e['by_label'], list), path + '.by_label')
        if e['retrieval_complete']:
            self.ok(not e['rows_refused'] and not e['anchors_truncated'], path,
                    'complete: nothing refused or truncated')
            self.ok(e['counts']['labels']['relation'] == 'exact', path + '.counts.labels')
            self.ok(e['counts']['occurrences']['relation'] == ('exact' if records else 'unknown'),
                    path + '.counts.occurrences')
        if not records:
            self.ok(e['counts']['occurrences']['relation'] == 'unknown',
                    path + '.counts.occurrences', 'placed occurrences need a record mapping')
        if e['counts']['labels']['relation'] == 'exact':
            self.ok(e['counts']['labels']['value'] == len(e['by_label'])
                    or e['labels_cut'] is not None, path + '.counts.labels',
                    'every label in by_label')
        keys = []
        for i, b in enumerate(e['by_label']):
            bp = f'{path}.by_label[{i}]'
            self.keys(b, SCHEMA['by_label'], bp)
            self.ok(b['graph'] == a['index']['index_ns'], bp + '.graph')
            self.count(b['contexts'], 'graph_contexts', bp + '.contexts')
            self.count(b['contexts_suffix'], 'graph_contexts', bp + '.contexts_suffix')
            self.count(b['occurrences'], 'placed_occurrences', bp + '.occurrences')
            self.ok(b['contexts']['value'] >= 1, bp + '.contexts')
            keys.append((-b['contexts']['value'], b['column']))
        self.ok(keys == sorted(keys), path + '.by_label', '(contexts desc, column asc)')
        rank = {b['column']: i for i, b in enumerate(e['by_label'])}
        for i, r in enumerate(e['results']):
            rp = f'{path}.results[{i}]'
            self.ok(r['support'] == 'kmer', rp + '.support')
            self.one_of(r['labels_status'], LABELS_STATUS, rp + '.labels_status')
            read = r['labels_status'] in ('complete', 'truncated')
            self.ok((r['labels_total'] is None) is (r['labels_status'] in ('refused', 'not_read')),
                    rp + '.labels_total')
            if not read:
                self.ok(r['labels'] is None, rp + '.labels', 'null unless read')
                continue
            self.ok(isinstance(r['labels'], list), rp + '.labels')
            listed = [label['column'] for label in r['labels']]
            self.ok(all(c in rank for c in listed), rp + '.labels', 'labels of by_label')
            self.ok(listed == sorted(listed, key=rank.get), rp + '.labels', 'label order')
            if r['labels_status'] == 'complete' and e['labels_cut'] is None:
                self.ok(len(listed) == r['labels_total'], rp + '.labels', 'every label')
            for j, label in enumerate(r['labels']):
                lp = f'{rp}.labels[{j}]'
                self.keys(label, ['column', 'support'] + (['occurrences'] if records else [])
                          + (['occurrence_list'] if placed else []), lp)
                self.ok(label['support'] == 'kmer', lp + '.support')
                if records:
                    self.count(label['occurrences'], 'placed_occurrences', lp + '.occurrences')
                if not placed or label['occurrence_list'] is None:
                    continue
                for n, o in enumerate(label['occurrence_list']):
                    op = f'{lp}.occurrence_list[{n}]'
                    strand = r.get('strand')
                    if records:
                        self.keys(o, SCHEMA['occurrence'], op)
                        m = re.fullmatch(r'(\d+)-(\d+)', o['nt_coords'])
                        self.ok(m and int(m.group(2)) == int(m.group(1)) + L - 1
                                and int(m.group(1)) >= 1 + r['offset']
                                and int(m.group(2)) <= o['nt_length'], op + '.nt_coords')
                    else:
                        self.keys(o, SCHEMA['occurrence_global'], op)
                        self.ok(o['offset'] == r['offset'], op + '.offset')
                    self.ok(o['strand'] == strand, op + '.strand', 'the context\'s')

    # ------------------------------------------------------------- capabilities

    def block(self, b, path, multi):
        if multi:
            self.keys(b, SCHEMA['capabilities_multi'], path)
            self.ok(b == {'pattern_contract_version': 1, 'available': False,
                          'unavailable_reason': 'multi_graph_later_increment'}, path)
            return
        self.keys(b, SCHEMA['capabilities'], path)
        self.ok(b['pattern_contract_version'] == 1, path + '.pattern_contract_version')
        self.ok(b['available'] in (True, False, None), path + '.available')
        if b['available'] is False:
            self.one_of(b['unavailable_reason'], UNAVAILABLE, path + '.unavailable_reason')
        else:
            self.ok(b['unavailable_reason'] is None, path + '.unavailable_reason')
        for f, values in (('modes', MODES), ('projections', ('none', 'all', 'predicate_only')),
                          ('kinds', ('dna', 'iupac', 'protein')),
                          ('strands', ('both', 'forward', 'reverse'))):
            self.ok(isinstance(b[f], list) and b[f] and set(b[f]) <= set(values), f'{path}.{f}')
        self.ok(b['modes'] == list(MODES), path + '.modes')
        self.ok('none' in b['projections'], path + '.projections', '"none" is served everywhere')
        self.ok(b['default_mode'] in b['modes'], path + '.default_mode')
        self.ok(b['default_projection'] in b['projections'], path + '.default_projection')
        self.ok(not set(b['projections']) & set(b['projections_later_increment']),
                path + '.projections_later_increment')
        self.ok(not set(b['kinds']) & set(b['kinds_later_increment']),
                path + '.kinds_later_increment')
        self.ok(b['default_scope'] in ('suffix', 'any_offset'), path + '.default_scope')
        self.ok(b['default_strands'] in b['strands'], path + '.default_strands')
        self.keys(b['scopes_by_graph_mode'], GRAPH_MODES, path + '.scopes_by_graph_mode')
        self.keys(b['caps'], CAPS, path + '.caps')
        self.ok(all(is_num(v) for v in b['caps'].values()), path + '.caps')
        self.ok(is_num(b['default_time_budget_ms'])
                and b['finalize_reserve_ms'] < b['default_time_budget_ms']
                <= b['caps']['time_budget_ms'], path + '.default_time_budget_ms')
        self.ok(isinstance(b['caps_rule'], str), path + '.caps_rule')
        self.ok(b['resident_only'] is True, path + '.resident_only')
        if b['available'] is None:
            for f in ('graph_mode', 'k', 'alphabet', 'strand_stated', 'mask', 'scopes',
                      'placement', 'support', 'annotation'):
                self.ok(b[f] is None, f'{path}.{f}', 'unknown while the index loads')
            return
        self.ok(is_int(b['k']), path + '.k')
        if b['graph_mode'] is None:
            return
        self.one_of(b['graph_mode'], GRAPH_MODES, path + '.graph_mode')
        self.ok(b['scopes'] == b['scopes_by_graph_mode'][b['graph_mode']], path + '.scopes')
        self.ok(b['strand_stated'] is (b['graph_mode'] == 'basic'), path + '.strand_stated')
        self.one_of(b['mask'], MASKS, path + '.mask')
        self.ok((b['mask'] == 'absent') is (b['unavailable_reason'] == 'mask_required'),
                path + '.mask', 'mask absent exactly when mask_required')
        for f, values in (('placement', PLACEMENTS), ('support', SUPPORTS),
                          ('annotation', ANNOTATIONS)):
            self.ok(b[f] is None or b[f] in values, f'{path}.{f}')
        if b['graph_mode'] != 'basic' and b['placement'] is not None:
            self.ok(b['placement'] == 'none_canonical', path + '.placement')


class TestPatternFixtures(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.index = load('index.json')
        cls.fixtures = cls.index['fixtures']
        cls.bodies = {name: (load(name, 'request.json'), load(name, 'answer.json'))
                      for name in cls.fixtures}

    def capabilities_of(self, server):
        """The pattern block of the server a fixture ran on (its GET /capabilities fixture)."""
        for name, f in self.fixtures.items():
            if f['server'] == server and f['method'] == 'GET':
                return self.bodies[name][1]['pattern']
        return None

    def test_spec_names_the_fields(self):
        if not os.path.isfile(SPEC):
            self.skipTest('docs/SPEC-pattern-search.md is not in this checkout')
        with open(SPEC, encoding='utf-8') as f:
            text = f.read()
        tables = spec_tables(text)
        self.assertEqual(set(SCHEMA), set(tables), 'the schema tables of the SPEC')
        for name, fields in SCHEMA.items():
            self.assertEqual(sorted(fields), sorted(tables[name]), f'SPEC table {name}')
            self.assertEqual(len(tables[name]), len(set(tables[name])), f'{name}: a field twice')
        # every value an answer may carry is named in the SPEC
        for values in (MODES, SCOPES, RELATIONS, UNITS, GRAPH_MODES, STOP_PHASES, STOP_REASONS,
                       WITHHELD, CUT, NOTES, SLOT_ERRORS, REFUSALS, UNAVAILABLE, MASKS,
                       PLACEMENTS, SUPPORTS, ANNOTATIONS, tuple(ABSENCE_SCOPES.values())):
            for v in values:
                self.assertIn(f'`{v}`', text, f'the SPEC does not name `{v}`')

    def test_files(self):
        listed = set(self.fixtures)
        on_disk = {d for d in os.listdir(DATA) if os.path.isdir(os.path.join(DATA, d))}
        self.assertEqual(listed, on_disk)
        self.assertEqual(1, self.index['pattern_contract_version'])
        with open(os.path.join(DATA, 'README.md'), encoding='utf-8') as f:
            readme = f.read()
        for name, f in self.fixtures.items():
            self.assertEqual(1, readme.count(f'- `{name}` '), f'README line of {name}')
            self.assertIn(f['server'], self.index['servers'])
            self.assertIn(f['method'], ('GET', 'POST'))
            request = self.bodies[name][0]
            if f['method'] == 'GET':
                self.assertIsNone(request, name)
            else:
                self.assertEqual('/pattern', f['path'], name)
                self.assertIsInstance(request, dict, name)

    def test_answers(self):
        for name, f in self.fixtures.items():
            request, answer = self.bodies[name]
            check = Checker(self, name)
            with self.subTest(fixture=name):
                if f['method'] == 'GET':
                    self.capabilities_document(check, f, answer)
                elif f['status'] == 200:
                    self.request_200(check, request)
                    check.answer(answer, request, self.capabilities_of(f['server']))
                else:
                    self.refusal(check, f, answer)

    def request_200(self, check, request):
        """A request the server answered uses only the SPEC's fields (and, but in mode count,
        output.labels "none" or "all")."""
        check.ok(set(request) <= set(SCHEMA['request']), 'request', sorted(request))
        for i, p in enumerate(request['patterns']):
            check.ok(set(p) <= set(SCHEMA['request_pattern']), f'request.patterns[{i}]')
        check.ok(set(request.get('output', {})) <= set(SCHEMA['request_output']),
                 'request.output')

    def refusal(self, check, f, answer):
        if f['status'] == 503 and 'code' not in answer:
            # the index is loading: every route's body, no code
            check.keys(answer, ['error'], 'answer')
            check.ok(f['headers'].get('Retry-After') == '60', 'headers')
            return
        check.keys(answer, SCHEMA['refusal'], 'answer')
        check.one_of(answer['code'], REFUSALS, 'answer.code')
        check.ok(isinstance(answer['error'], str) and answer['error'], 'answer.error')
        check.ok((f['status'] == 503) is (answer['code'] == 'deadline'), 'status',
                 '503 is the deadline, every other refusal 400')
        check.ok(f['status'] in (400, 503), 'status')

    def capabilities_document(self, check, f, doc):
        multi = f['server'] == 'multi'
        check.block(doc['pattern'], 'pattern', multi)
        if f['path'] == '/capabilities':
            check.ok(doc['mode'] == ('multi' if multi else 'single'), 'mode')
            listed = not multi
            check.ok(('pattern' in doc['features']) is listed, 'features')
            check.ok(doc['routes'].get('pattern') == ('POST /pattern' if listed else None),
                     'routes.pattern')
        else:
            check.ok(f['path'].startswith('/traverse/capabilities'), 'path')

    def test_the_block_is_the_same_on_both_routes(self):
        by_server = {}
        for name, f in self.fixtures.items():
            if f['method'] == 'GET':
                by_server.setdefault(f['server'], []).append(self.bodies[name][1]['pattern'])
        for server, blocks in by_server.items():
            for b in blocks[1:]:
                self.assertEqual(blocks[0], b, server)

    def test_requests_and_answers_agree(self):
        """What a fixture's request asks is what its answer states."""
        for name, f in self.fixtures.items():
            request, answer = self.bodies[name]
            if f['method'] != 'POST' or f['status'] != 200:
                continue
            with self.subTest(fixture=name):
                for e in answer['patterns']:
                    if 'error' in e or e['scope'] == 'long':
                        continue
                    if request.get('stop_at_threshold') and e['stop'] \
                            and e['stop']['reason'] == 'max_contexts':
                        c = e['counts']['contexts']
                        self.assertEqual('at_least', c['relation'])
                        self.assertGreater(c['value'], answer['limits']['max_contexts'])
                    if e.get('withheld') == {'reason': 'count_above_threshold'}:
                        c = e['counts']['contexts']
                        self.assertEqual('exact', c['relation'])
                        self.assertGreater(c['value'], answer['limits']['max_contexts'])
                    if e.get('cut') == {'reason': 'max_contexts'}:
                        self.assertEqual(answer['limits']['max_contexts'], e['returned'])

    def test_coverage(self):
        """The fixtures show every situation a client of version 1 has to handle."""
        seen = {k: set() for k in ('mode', 'scope', 'strands', 'withheld', 'cut', 'stop',
                                   'slot', 'refusal', 'graph_mode', 'relation', 'note')}
        for name, f in self.fixtures.items():
            request, answer = self.bodies[name]
            if f['method'] != 'POST':
                continue
            if f['status'] != 200:
                seen['refusal'].add(answer.get('code', 'initializing'))
                continue
            seen['mode'].add(answer['mode'])
            seen['strands'].add(request.get('strands', 'both'))
            seen['graph_mode'].add(answer['index']['graph_mode'])
            for e in answer['patterns']:
                if 'error' in e:
                    seen['slot'].add(e['error']['code'])
                    continue
                seen['scope'].add(e['scope'])
                seen['note'].update(e['notes'])
                for f2 in ('withheld', 'cut'):
                    if e.get(f2):
                        seen[f2].add(e[f2]['reason'])
                if e['stop']:
                    seen['stop'].add((e['stop']['phase'], e['stop']['reason']))
                total = e['counts'].get('contexts') or e['counts']['anchors']
                seen['relation'].add(total['relation'])
        self.assertEqual(set(MODES), seen['mode'])
        self.assertEqual(set(SCOPES), seen['scope'])
        self.assertEqual({'both', 'forward', 'reverse'}, seen['strands'])
        # output_budget and the cut max_memory need more contexts than the mini's patterns
        # give at max_memory_mb 1 (the unit tests show them, tests/cli/test_pattern_retrieval)
        self.assertEqual(set(WITHHELD) - {'output_budget'}, seen['withheld'])
        self.assertEqual({'max_contexts', 'max_steps', 'time'}, seen['cut'])
        self.assertEqual(set(SLOT_ERRORS), seen['slot'])
        self.assertEqual({'basic', 'primary'}, seen['graph_mode'])
        self.assertEqual(set(NOTES), seen['note'])
        self.assertLessEqual({'exact', 'at_least', 'unknown'}, seen['relation'])
        self.assertLessEqual({('discovery', 'max_steps'), ('discovery', 'time'),
                              ('discovery', 'max_contexts')}, seen['stop'])
        self.assertEqual({'invalid_request', 'later_increment', 'resident_only', 'mask_required',
                          'deadline', 'initializing', 'annotation_unbudgeted'}, seen['refusal'])
        self.assertLessEqual({('label_discovery', 'max_annotation_work')}, seen['stop'])

    def test_hand_made_bodies_are_the_codes(self):
        """The two hand-made 503 bodies are written by the code as stored here."""
        if not (os.path.isfile(PATTERN_CPP) and os.path.isfile(SERVER_UTILS_CPP)):
            self.skipTest('the server sources are not in this checkout')
        with open(PATTERN_CPP, encoding='utf-8') as f:
            pattern_cpp = f.read()
        with open(SERVER_UTILS_CPP, encoding='utf-8') as f:
            server_utils = f.read()
        hand_made = {n for n, f in self.fixtures.items() if f['hand_made']}
        self.assertEqual({'deadline_503', 'initializing_503'}, hand_made)

        _, deadline = self.bodies['deadline_503']
        self.assertEqual('deadline', deadline['code'])
        m = re.fullmatch(r'pattern: the answer could not be written within time_budget_ms '
                         r'\((\d+) ms, the finalisation reserve of (\d+) ms included\): '
                         r'nothing partial is sent', deadline['error'])
        self.assertTrue(m, deadline['error'])
        for piece in ('"pattern: the answer could not be written within time_budget_ms ("',
                      '" ms, the finalisation reserve of "', '" ms included): nothing partial '
                      'is sent"', 'PatternRefusal(503, "deadline"'):
            self.assertIn(piece, pattern_cpp)
        self.assertEqual(self.bodies['deadline_503'][0]['time_budget_ms'], int(m.group(1)))

        _, loading = self.bodies['initializing_503']
        self.assertIn(json.dumps(loading['error']), server_utils)
        self.assertIn('"Retry-After", "60"', server_utils)

    def test_mask_required_text_is_the_codes(self):
        """The mask_required refusal stored is the text pattern.cpp writes now (review of
        milestone 1b: the fixtures had been generated by the milestone-1 binary, whose message
        named only `build --mask-dummy`; this check needs no server)."""
        if not os.path.isfile(PATTERN_CPP):
            self.skipTest('the server sources are not in this checkout')
        with open(PATTERN_CPP, encoding='utf-8') as f:
            pattern_cpp = f.read()
        # the C++ string literals of support_message's mask_required branch, concatenated
        m = re.search(r'support\.reason == "mask_required"\) \{\s*return ((?:"(?:[^"\\]|\\.)*"'
                      r'\s*)+);', pattern_cpp)
        self.assertTrue(m, 'support_message\'s mask_required branch not found in pattern.cpp')
        text = ''.join(json.loads(lit) for lit in re.findall(r'"(?:[^"\\]|\\.)*"', m.group(1)))
        _, answer = self.bodies['mask_required']
        self.assertEqual({'error': text, 'code': 'mask_required'}, answer)
        # both remedies of DESIGN §4: the file once, or the mask built at every start
        self.assertIn('metagraph transform --mask-dummy', text)
        self.assertIn('--pattern-build-mask', text)

    def test_every_mask_value_has_a_capabilities_fixture(self):
        """file, built_at_load and absent (DESIGN §4), each on both capabilities routes."""
        seen = {}
        for name, f in self.fixtures.items():
            if f['method'] == 'GET' and f['server'] != 'multi':
                route = 'probe' if f['path'].startswith('/traverse/') else 'capabilities'
                seen.setdefault(self.bodies[name][1]['pattern']['mask'], set()).add(route)
        self.assertEqual(set(MASKS), set(seen))
        for mask in ('file', 'built_at_load', 'absent'):
            self.assertEqual({'probe', 'capabilities'}, seen[mask], mask)

    def test_documents_state_the_built_deadline(self):
        """The route's default deadline and its cap, as the SPEC and the service's request
        (PROMPT-search-service-pattern.md §1 and §3.2) state them, are the capabilities'
        (review of milestone 1b: PROMPT §1 still said "default 5 s")."""
        caps = self.bodies['capabilities'][1]['pattern']
        default_s = caps['default_time_budget_ms'] / 1000
        cap_s = caps['caps']['time_budget_ms'] / 1000
        self.assertEqual((60, 600), (default_s, cap_s))
        if os.path.isfile(PROMPT):
            with open(PROMPT, encoding='utf-8') as f:
                prompt = f.read()
            stated = re.findall(r'its own deadline \(default (\d+) s, cap\s+(\d+) s', prompt)
            self.assertEqual([(str(int(default_s)), str(int(cap_s)))], stated, 'PROMPT §1')
            self.assertEqual([str(int(default_s))],
                             re.findall(r'\*\*default\*\* `time_budget_ms` is \*\*(\d+) s\*\*',
                                        prompt), 'PROMPT §3.2')
            self.assertEqual([str(int(cap_s))],
                             re.findall(r'`--pattern-max-time-ms` is \*\*(\d+) s\*\*', prompt),
                             'PROMPT §3.2')
            self.assertNotRegex(prompt, r'default\s+5\s+s\b')
        if os.path.isfile(SPEC):
            with open(SPEC, encoding='utf-8') as f:
                spec = f.read()
            for flag, value in (('--pattern-default-time-ms', caps['default_time_budget_ms']),
                                ('--pattern-max-time-ms', caps['caps']['time_budget_ms'])):
                self.assertIn(f'| `{flag}` | {value:,} |', spec, flag)


if __name__ == '__main__':
    unittest.main()
