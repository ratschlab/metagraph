"""The POST /pattern fixture bodies (data/traverse/pattern/, written by
scripts/traversal/pattern_fixtures.py) against docs/SPEC-pattern-search.md, contract version 1
(milestone 1, increment 3's output.labels "all", SPEC §14, and increments 4 and 5, SPEC §12.1,
§12.2, §17: the paths of long_search "paths" with their labels' support, and protein patterns):

  - the field lists: every table of the SPEC marked `<!-- schema: NAME -->` names exactly the
    fields SCHEMA[NAME] below knows, so that the SPEC and this check cannot drift apart;
  - every answer: exactly the fields its shape has (answered entry, error slot, refusal,
    capabilities block), the types and enumerations of the SPEC, and the rules that tie the
    fields together -- a count's relation and its value, a total and its parts (exact: their
    sum; at_least and bounds: the sum of their lower bounds, bounds also of their upper ones;
    the weakest relation), the offsets of the scope, the information floor and its exemption,
    a clamp to the cap, one budget per request (the steps of all patterns within max_steps, a
    max_steps or time stop sticky: every later pattern unknown, work zero, the same stop), the
    order of the results and, for a complete list, its agreement with the counts per offset and
    per strand, what withheld, cut and retrieval_complete imply -- checked on the bodies, never
    by asking a server (review of 2026-10-07, C2-02, C2-03, X-ORACLE-03: the cross-field rules
    were missing); a flank of a result's k-mer is checked against the graph by
    `pattern_fixtures.py --check` only;
  - the fixture set: it covers every mode, scope, strand setting, relation (bounds included),
    stop, withheld and cut reason, slot error, refusal code and unavailable reason a service
    has to handle -- but for the ones UNPRODUCIBLE names, which no server of this build can
    give, and the ones NO_FIXTURE names, which need an index the fixture servers are not --
    and the README names each fixture; the refusal codes and unavailable reasons are the ones
    the server's sources write.

No metagraph import and no server: the service copies these bodies into its own tests, and
this check is what makes them a contract rather than a sample. A peptide's bits and instances
are checked against the validator's own copy of NCBI's genetic codes (GENETIC_CODES), never
against the server's tables.
"""

import copy
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
# increment 3 (output.labels "all", SPEC §14) adds the annotation phases and reasons;
# increment 4 (long_search "paths", SPEC §17) the phase extension and the reason max_paths
ANNOTATION_PHASES = ('label_discovery', 'placement', 'output')
STOP_PHASES = ('discovery', 'mask_scan', 'extraction') + ANNOTATION_PHASES + ('extension',)
STOP_REASONS = ('max_steps', 'time', 'max_contexts', 'max_anchors', 'max_annotation_work',
                'max_memory', 'max_paths')
WITHHELD = ('count_above_threshold', 'threshold_crossed', 'discovery_budget', 'deadline',
            'paths_later_increment', 'annotation_budget', 'anchor_labels_truncated',
            'output_budget', 'anchors_above_threshold')
CUT = ('max_contexts', 'max_steps', 'time', 'max_anchors', 'max_memory', 'max_paths')
NOTES = ('low_complexity_pattern', 'strand_unknown_canonical', 'paths_later_increment',
         'annotation_unbudgeted', 'record_bounds_unknown', 'annotation_not_read',
         'label_intersection_only')
# stop_unsupported: a peptide valid apart from its stop '*' (increment 5, owner decision #15)
SLOT_ERRORS = ('bad_alphabet', 'information_below_floor', 'scope_unsupported',
               'stop_unsupported')
KINDS = ('dna', 'iupac', 'protein')
# increment 4: what the extension of a pattern longer than k did (counts.paths.extension)
EXTENSIONS = ('no_anchors', 'not_started', 'not_admitted', 'stopped', 'completed')
LONG_SEARCH = ('anchors', 'paths')
# a path result's support: its listed labels' (null when it lists none)
PATH_SUPPORTS = ('record_verified', 'label_intersection', 'mixed')
# increment 5: the residues of a protein pattern (any case) and NCBI's genetic codes (gc.prt
# version 4.6, https://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt): each table's ncbieaa string,
# the residue of each codon in TCAG order (TTT TTC TTA TTG TCT ... GGG), '*' a stop
PROTEIN_RESIDUES = 'ACDEFGHIKLMNPQRSTVWYXBZJ'
GENETIC_CODES = {
    1: 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    2: 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG',
    3: 'FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    4: 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    5: 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG',
    6: 'FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    9: 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG',
    10: 'FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    11: 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    12: 'FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    13: 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG',
    14: 'FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG',
    15: 'FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    16: 'FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    21: 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG',
    22: 'FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    23: 'FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    24: 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG',
    25: 'FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    26: 'FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    27: 'FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    28: 'FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    29: 'FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    30: 'FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    31: 'FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    32: 'FFLLSSSSYY*WCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    33: 'FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG',
}
CODONS_TCAG = [a + b + c for a in 'TCAG' for b in 'TCAG' for c in 'TCAG']
# mask_invalid and alphabet_untested: the route's own reasons not to serve a graph (85614d30,
# owner decisions #6 and #4 of 2026-10-07; pattern.cpp route_support), a 400 refusal's code and
# the capabilities' unavailable_reason alike (review GPT-2 of 2026-10-08, finding 6: missing
# here, valid answers carrying them were rejected)
REFUSALS = ('invalid_request', 'later_increment', 'resident_only', 'mask_required',
            'representation_unsupported', 'primary_unwrapped', 'alphabet_unsupported', 'deadline',
            'annotation_unbudgeted', 'mask_invalid', 'alphabet_untested',
            # increment 4 (require_support "record_verified" on an index that cannot verify) and
            # increment 5 (genetic_code no NCBI table)
            'support_unavailable', 'genetic_code_unknown')
UNAVAILABLE = ('mask_required', 'representation_unsupported', 'primary_unwrapped',
               'alphabet_unsupported', 'multi_graph_later_increment', 'mask_invalid',
               'alphabet_untested')
# The refusals and unavailable reasons no server of this build gives, so that no fixture holds
# them (review of 2026-10-07, C2-01, D1-07, C1-05): primary_unwrapped (server_query and the CLI
# always wrap a PRIMARY graph in CanonicalDBG; only an embedding reaches it) and
# alphabet_unsupported (the BOSS alphabet is the build's: a graph of another alphabet does not
# load).
UNPRODUCIBLE = ('primary_unwrapped', 'alphabet_unsupported')
# The ones a server of this build gives but no stored fixture holds, each with the index it
# needs (review GPT-2 of 2026-10-08, finding 6): named here rather than left out of the
# expected sets, so that the coverage tests list every code without a fixture. Every other code
# has a fixture.
NO_FIXTURE = {
    'mask_invalid': 'a graph whose .edgemask marks an edge with W = $ valid: one extended after '
                    'masking by a build older than 85614d30 (metagraph extend on a masked '
                    'graph), or a stale mask',
    'alphabet_untested': 'a DNA5 build of metagraph and a $ACGTN graph (a DNA4 build does not '
                         'load one)',
}
# The alphabet the route serves (pattern.cpp alphabet_refusal): on any other the alphabet's
# reason comes first, mask or none
SERVED_ALPHABET = '$ACGT'
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
                'max_occurrences_per_label', 'allow_unbudgeted_annotation',
                'long_search', 'max_paths', 'require_support', 'genetic_code'],
    'request_pattern': ['id', 'dna', 'iupac', 'protein'],
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
               'allow_unbudgeted_annotation', 'long_search', 'max_paths', 'require_support',
               'clamped'],
    'clamped': ['field', 'requested', 'effective'],
    'timing': ['elapsed_ms', 'label_discovery_ms', 'placement_ms', 'extension_ms'],
    'entry': ['id', 'kind', 'pattern', 'length', 'residues', 'genetic_code', 'information_bits',
              'anchor_information_bits', 'min_anchor_information_bits', 'error', 'mode', 'scope',
              'strands', 'palindromic', 'counts', 'work', 'stop',
              'retrieval_complete', 'withheld', 'returned', 'cut', 'results', 'absence_scope',
              'determinism', 'notes', 'timing', 'placement', 'annotation', 'by_label',
              'rows_refused', 'anchors_truncated', 'labels_cut', 'occurrences_cut',
              'labels_excluded_unverified'],
    'counts': ['contexts', 'anchors', 'paths', 'labels', 'occurrences'],
    'count': ['value', 'relation', 'unit', 'lower', 'upper'],
    'contexts_count': ['suffix', 'by_offset', 'by_strand', 'by_orientation'],
    'anchors_count': ['by_strand', 'by_orientation'],
    # increment 4 (long_search "paths", SPEC §17)
    'paths_count': ['by_strand', 'by_orientation', 'candidates_examined', 'extension'],
    'labels_count': ['by_support'],
    'by_support': ['record_verified', 'label_intersection'],
    'work': ['ranges_visited', 'mask_scans', 'steps', 'annotation_rows', 'annotation_units',
             'memory_bytes', 'extension_edges'],
    'stop': ['phase', 'reason'],
    'reason': ['reason'],
    'result': ['kmer', 'instance', 'offset', 'strand', 'orientation', 'node', 'row', 'support',
               'labels_status', 'labels_total', 'labels'],
    # increment 4: a path (a result of a pattern longer than k under long_search "paths")
    'path_result': ['sequence', 'anchor_kmer', 'instance', 'offset', 'strand', 'orientation',
                    'nodes', 'rows', 'support', 'labels_status', 'labels_total', 'labels',
                    'labels_excluded_unverified'],
    'error': ['code', 'message'],
    'capabilities': ['pattern_contract_version', 'available', 'unavailable_reason', 'modes',
                     'default_mode', 'projections', 'default_projection',
                     'projections_later_increment', 'kinds', 'kinds_later_increment',
                     'default_scope', 'scopes_by_graph_mode', 'scopes', 'long_patterns',
                     'strands', 'default_strands', 'graph_cleaned', 'records_shorter_than_k',
                     'resident_only', 'caps', 'default_time_budget_ms', 'finalize_reserve_ms',
                     'caps_rule', 'graph_mode', 'k', 'alphabet', 'strand_stated', 'mask',
                     'placement', 'support', 'annotation', 'default_occurrences',
                     # increments 4 and 5
                     'long_search', 'default_long_search', 'protein_residues', 'genetic_codes',
                     'default_genetic_code', 'protein_rule'],
    'capabilities_multi': ['pattern_contract_version', 'available', 'unavailable_reason'],
    # increment 3 (SPEC §14)
    'anchor_truncated': ['kmer', 'row', 'cap', 'total'],
    'occurrence': ['seq_id', 'record', 'strand', 'nt_coords', 'nt_length'],
    'occurrence_global': ['kmer_coord', 'offset', 'strand'],
    'row_refused': ['kmer', 'row', 'phase', 'reason', 'needed_bytes', 'available_bytes'],
    'label': ['column', 'support', 'occurrences', 'occurrence_list'],
    'by_label': ['graph', 'column', 'contexts', 'contexts_suffix', 'occurrences'],
    'by_label_paths': ['graph', 'column', 'paths', 'paths_record_verified', 'occurrences'],
    'labels_cut': ['reason', 'returned'],
    'occurrences_cut': ['reason', 'labels'],
}
ENTRY_DESCRIPTION = ['id', 'kind', 'pattern', 'length', 'information_bits',
                     'anchor_information_bits', 'min_anchor_information_bits']
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
        'max_patterns'] + LIMITS_LABELS + ['max_paths']


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


def codons_of(residue, table):
    """The codons of a peptide's residue in an NCBI table (X: every codon that is not a stop,
    B: D or N, Z: E or Q, J: I or L; no stop codon ever)."""
    admits = {'X': lambda a: True, 'B': lambda a: a in 'DN', 'Z': lambda a: a in 'EQ',
              'J': lambda a: a in 'IL'}.get(residue, lambda a: a == residue)
    return {c for c, a in zip(CODONS_TCAG, GENETIC_CODES[table]) if a != '*' and admits(a)}


class Kind:
    """A request pattern as the validator reads it: dna or iupac position by position, a
    protein (a peptide) as its codon sets in the request's genetic code (SPEC §12.2)."""

    def __init__(self, asked, request):
        self.kind = next(k for k in KINDS if k in asked)
        self.text = asked[self.kind].upper()
        self.table = request.get('genetic_code', 1)
        self.protein = self.kind == 'protein'
        if self.protein and self.parsed():
            self.codons = [codons_of(r, self.table) for r in self.text]

    def parsed(self):
        """Within the kind's alphabet (no slot error of the alphabet)."""
        if not self.text:
            return False
        alphabet = {'dna': 'ACGT', 'iupac': ''.join(IUPAC), 'protein': PROTEIN_RESIDUES}
        return all(c in alphabet[self.kind] for c in self.text)

    def length(self):
        return 3 * len(self.text) if self.protein else len(self.text)

    def bits(self, begin=0, end=None):
        """The information of the bases [begin, end): DNA and IUPAC per position; a peptide per
        residue, 2 per base less log2 of the distinct strings its codons spell over the bases
        of the window it covers (exact also for a window that cuts a codon, SPEC §12.2)."""
        end = self.length() if end is None else end
        if not self.protein:
            return information_bits(self.text[begin:end])
        out = 0.0
        for i, codons in enumerate(self.codons):
            lo, hi = max(begin, 3 * i) - 3 * i, min(end, 3 * i + 3) - 3 * i
            if lo < hi:
                out += 2 * (hi - lo) - math.log2(len({c[lo:hi] for c in codons}))
        return out

    def instance(self, s, reverse):
        """|s| is an instance of the oriented pattern: P, or rc(P) when |reverse|."""
        if not self.protein:
            return matches(revcomp(self.text) if reverse else self.text, s)
        m = len(self.text)
        if len(s) != 3 * m:
            return False
        for i in range(m):
            codon = s[3 * i:3 * i + 3]
            if reverse:
                if revcomp(codon) not in self.codons[m - 1 - i]:
                    return False
            elif codon not in self.codons[i]:
                return False
        return True

    def palindromic(self):
        if not self.protein:
            return revcomp(self.text) == self.text
        m = len(self.codons)
        return all(self.codons[i] == {revcomp(c) for c in self.codons[m - 1 - i]}
                   for i in range(m))

    def exact(self):
        """Every position one base (a peptide: every residue one codon)."""
        if self.protein:
            return all(len(c) == 1 for c in self.codons)
        return all(c in 'ACGT' for c in self.text)


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
        # SPEC §7.4 (Count::operator+=): an at_least or bounds total is the sum of the lower
        # bounds of what was explored (unknown as 0), a bounds total's upper the sum of the
        # uppers; the relation is the weakest of the parts'
        def low(c):
            return 0 if c['relation'] == 'unknown' else c['value']

        def high(c):
            return c['upper'] if c['relation'] == 'bounds' else c['value']
        if total['relation'] in ('at_least', 'bounds'):
            self.ok(total['value'] == sum(low(p) for p in parts), path,
                    'an at_least or bounds total is the sum of its parts\' lower bounds')
            self.ok(any(p['relation'] in ('at_least', 'unknown') for p in parts)
                    is (total['relation'] == 'at_least'), path, 'the weakest relation')
        if total['relation'] == 'bounds':
            self.ok(total['upper'] == sum(high(p) for p in parts), path,
                    'a bounds total\'s upper is the sum of its parts\' uppers')

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
        # the annotation limits only in the answers that read annotation (increment 3); the
        # path limits only in the answers to long_search "paths" (increment 4, SPEC §12.1),
        # require_support among them only when labels are read
        paths = request.get('long_search', 'anchors') == 'paths'
        self.keys(limits, [f for f in SCHEMA['limits']
                           if (labelled or f not in LIMITS_LABELS + ['allow_unbudgeted_annotation'])
                           and (paths or f not in ('long_search', 'max_paths'))
                           and (paths and labelled or f != 'require_support')],
                  'limits')
        if labelled:
            for f in LIMITS_LABELS:
                self.ok(is_int(limits[f]) and limits[f] >= 0, 'limits.' + f)
            self.ok(limits['allow_unbudgeted_annotation']
                    is request.get('allow_unbudgeted_annotation', False),
                    'limits.allow_unbudgeted_annotation')
        if paths:
            self.ok(limits['long_search'] == 'paths', 'limits.long_search')
            self.ok(is_int(limits['max_paths']) and limits['max_paths'] >= 0, 'limits.max_paths')
            if labelled:
                self.ok(limits['require_support']
                        == request.get('require_support', 'label_intersection'),
                        'limits.require_support', 'as requested')
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
        order = ['max_contexts', 'max_anchors', 'max_steps', 'time_budget_ms'] + LIMITS_LABELS \
            + ['max_paths']
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
            if capabilities is not None:
                # SPEC §4.5: lowered to the cap, not to the default
                self.ok(c['effective'] == capabilities['caps'][c['field']], path,
                        'lowered to the cap')
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
            if 'genetic_code' in request:
                self.ok(request['genetic_code'] in capabilities['genetic_codes'],
                        'request.genetic_code', 'a table of the capabilities')

        self.timing(a['timing'], 'timing')
        self.ok(len(a['patterns']) == len(request['patterns']), 'patterns',
                'one entry per request pattern')
        for i, (entry, asked) in enumerate(zip(a['patterns'], request['patterns'])):
            self.entry(entry, asked, request, a, f'patterns[{i}]')
        self.budget(a, limits)

    def budget(self, a, limits):
        """SPEC §7.6, one budget per request: the steps of all patterns sum to at most max_steps,
        and a stop by max_steps or time is sticky -- every later answered pattern stops with the
        same reason in discovery, its counts unknown, its (engine) work zero."""
        answered = [e for e in a['patterns'] if 'error' not in e]
        self.ok(sum(e['work']['steps'] for e in answered) <= limits['max_steps'], 'patterns',
                'the steps of all patterns sum to at most max_steps')
        sticky = None
        timed = False
        for i, e in enumerate(a['patterns']):
            if 'error' in e:
                continue
            path = f'patterns[{i}]'
            # SPEC §7.9: the entry a time stop touched and every answered entry after it
            if timed:
                self.ok(e['determinism'] == 'time_limited', path + '.determinism',
                        'an entry after a time-limited one is time-limited')
            timed = timed or e['determinism'] == 'time_limited'
            if sticky is not None:
                self.ok(e['stop'] == {'phase': 'discovery', 'reason': sticky}, path + '.stop',
                        f'the {sticky} stop of an earlier pattern, in discovery')
                engine = ('ranges_visited', 'mask_scans', 'steps', 'extension_edges')
                self.ok(all(e['work'][f] == 0 for f in engine if f in e['work']), path + '.work',
                        'no work after a budget stop')
                for name, c in self.counts_of(e):
                    self.ok(c['relation'] == 'unknown' and c['value'] is None,
                            f'{path}.counts.{name}', 'unknown after a budget stop')
            elif e['stop'] is not None and e['stop']['reason'] in ('max_steps', 'time'):
                sticky = e['stop']['reason']

    @staticmethod
    def counts_of(e):
        """Every graph count of an answered entry: (name, count), parts included."""
        c = e['counts']
        out = []
        for unit in ('contexts', 'anchors', 'paths'):
            if unit not in c:
                continue
            top = c[unit]
            out.append((unit, top))
            if 'suffix' in top:
                out.append((unit + '.suffix', top['suffix']))
            for key in ('by_offset', 'by_strand', 'by_orientation'):
                for k, v in top.get(key, {}).items():
                    out.append((f'{unit}.{key}.{k}', v))
        return out

    def timing(self, t, path, labelled=False, extension=False):
        self.keys(t, ['elapsed_ms'] + (['label_discovery_ms', 'placement_ms'] if labelled else [])
                  + (['extension_ms'] if extension else []), path)
        for k, v in t.items():
            self.ok(is_num(v) and v >= 0, f'{path}.{k}')

    def paths_count(self, c, anchors, e, strand_stated, limits, path):
        """counts.paths of a pattern longer than k under long_search "paths" (SPEC §12.1): the
        count with its per-orientation split, candidates_examined and what the extension did,
        each extension with its relation."""
        key = 'by_strand' if strand_stated else 'by_orientation'
        self.count(c, 'paths', path, extra=[key, 'candidates_examined', 'extension'])
        self.by_orientation(c, e, strand_stated, 'paths', path)
        self.ok(is_int(c['candidates_examined']) and c['candidates_examined'] >= 0,
                path + '.candidates_examined')
        self.one_of(c['extension'], EXTENSIONS, path + '.extension')
        relation = {'completed': 'exact', 'no_anchors': 'exact', 'stopped': 'at_least',
                    'not_started': 'unknown', 'not_admitted': 'unknown'}.get(c['extension'])
        self.ok(c['relation'] == relation, path + '.relation',
                f'{c["relation"]} for an extension {c["extension"]}')
        no_anchor = (anchors['relation'], anchors['value']) == ('exact', 0)
        self.ok((c['extension'] == 'no_anchors') is no_anchor, path + '.extension',
                'no_anchors exactly when the anchors are exact 0')
        if no_anchor:
            self.ok(c['value'] == 0 and c['candidates_examined'] == 0, path,
                    'no anchor, no path: exact 0')
        if c['extension'] == 'completed':
            self.ok(anchors['relation'] == 'exact', path + '.extension',
                    'completed: every anchor discovered and extended')
        if c['extension'] == 'not_admitted':
            self.ok(anchors['relation'] == 'exact' and anchors['value'] > limits['max_anchors'],
                    path + '.extension', 'not admitted: the exact anchors above max_anchors')
        if c['relation'] in ('exact', 'at_least'):
            # each path entered its last branch of its own
            self.ok(c['value'] <= c['candidates_examined'], path + '.candidates_examined',
                    'at least one branch per path')

    def by_support(self, labels, path):
        """counts.labels.by_support of the labels of paths: the labels verified on a returned
        path, and the listed ones verified on none; the two sum to the labels when exact."""
        b = labels['by_support']
        self.keys(b, SCHEMA['by_support'], path + '.by_support')
        for name, v in b.items():
            self.count(v, 'labels', f'{path}.by_support.{name}')
        if labels['relation'] == 'exact' and all(v['relation'] == 'exact' for v in b.values()):
            self.ok(sum(v['value'] for v in b.values()) == labels['value'], path + '.by_support',
                    'the supports sum to the labels')
        if labels['relation'] == 'unknown':
            self.ok(all(v['relation'] == 'unknown' for v in b.values()), path + '.by_support',
                    'unknown labels, unknown supports')

    def entry(self, e, asked, request, a, path):
        k = a['index']['k']
        strand_stated = a['index']['strand_stated']
        mode = a['mode']
        labelled = isinstance(a['output'], dict) and a['output']['labels'] == 'all'
        # the request named the labels (or an annotation field) and this answer reads none
        named = (request.get('output', {}).get('labels') == 'all'
                 or any(f in request for f in LIMITS_LABELS + ['allow_unbudgeted_annotation',
                                                               'require_support']))
        self.ok(set(e) <= set(SCHEMA['entry']), path, f'fields outside the SPEC: '
                f'{sorted(set(e) - set(SCHEMA["entry"]))}')
        self.ok(e['id'] == asked.get('id'), path + '.id', 'echoed')
        pk = Kind(asked, request)
        self.ok(e['kind'] == pk.kind, path + '.kind')
        text = pk.text
        # a peptide's entry states its residues and its genetic code (increment 5)
        description = ENTRY_DESCRIPTION + (['residues', 'genetic_code'] if pk.protein else [])

        if 'error' in e:
            err = e['error']
            self.keys(err, SCHEMA['error'], path + '.error')
            self.one_of(err['code'], SLOT_ERRORS, path + '.error.code')
            self.ok(isinstance(err['message'], str) and err['message'], path + '.error.message')
            # a peptide valid apart from its stops '*' (SPEC §8.9, §12.2)
            stops_only = pk.protein and '*' in text \
                and all(c in PROTEIN_RESIDUES + '*' for c in text)
            if err['code'] in ('bad_alphabet', 'stop_unsupported'):
                self.keys(e, ['id', 'kind', 'error'], path)
                if err['code'] == 'stop_unsupported':
                    self.ok(stops_only, path, 'stop_unsupported: a peptide valid but its stops')
                else:
                    self.ok(not pk.parsed() and not stops_only, path,
                            'bad_alphabet for a pattern inside the alphabet')
                return
            self.ok(pk.parsed(), path, 'an engine slot for a pattern outside the alphabet')
            self.keys(e, description + ['error'], path)
            self.description(e, pk, k, path, request)
            if err['code'] == 'information_below_floor':
                bits = e['min_anchor_information_bits'] if e['length'] > k \
                    else e['information_bits']
                self.ok(bits < a['limits']['min_information_bits'], path,
                        'refused above the floor')
                # SPEC §7.8: an exact pattern in suffix scope is exempt from the floor
                self.ok(not self.exempt(pk, k, request), path,
                        'an exact pattern of L <= k in suffix scope is exempt from the floor')
            return

        self.ok(pk.parsed(), path, 'answered outside the alphabet')
        L = pk.length()
        # long_search "paths" (increment 4): a pattern longer than k answered by its paths
        paths = request.get('long_search', 'anchors') == 'paths' and L > k
        verified_only = request.get('require_support') == 'record_verified'
        expected = description + ENTRY_ANSWERED[len(ENTRY_DESCRIPTION):] \
            + (ENTRY_RETRIEVAL if mode != 'count' else []) + (ENTRY_LABELS if labelled else []) \
            + (['labels_excluded_unverified'] if labelled and paths and verified_only else [])
        self.keys(e, expected, path)
        self.description(e, pk, k, path, request)
        # SPEC §7.8: a pattern answered is exempt or at or above the floor (for L > k: every
        # searched anchor window, review of 2026-10-07, X-GUARANTEES-01)
        bits = e['min_anchor_information_bits'] if L > k else e['information_bits']
        self.ok(self.exempt(pk, k, request) or bits >= a['limits']['min_information_bits'],
                path, 'answered below the floor')
        self.ok(e['mode'] == mode, path + '.mode')
        scope = 'long' if L > k else request.get('scope', 'any_offset')
        self.ok(e['scope'] == scope, path + '.scope', f'{e["scope"]!r}, expected {scope!r}')
        self.ok(e['absence_scope'] == ABSENCE_SCOPES[scope], path + '.absence_scope')
        palindromic = pk.palindromic()
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
            if paths:
                self.paths_count(counts['paths'], anchors, e, strand_stated, a['limits'],
                                 path + '.counts.paths')
            else:
                self.count(counts['paths'], 'paths', path + '.counts.paths')
                if (anchors['relation'], anchors['value']) == ('exact', 0):
                    self.ok(counts['paths'] == {'value': 0, 'relation': 'exact', 'unit': 'paths'},
                            path + '.counts.paths', 'no anchor, no path: exact 0')
                else:
                    self.ok(counts['paths']['relation'] == 'unknown', path + '.counts.paths',
                            'paths only with long_search "paths"')
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
            if labelled and f == 'labels' and paths:
                # the labels of paths: split by their support (SPEC §12.1)
                self.count(counts[f], unit, f'{path}.counts.{f}', extra=['by_support'])
                self.by_support(counts[f], f'{path}.counts.{f}')
            elif labelled:
                self.count(counts[f], unit, f'{path}.counts.{f}')
            else:
                self.ok(counts[f] == {'value': None, 'relation': 'unknown', 'unit': unit},
                        f'{path}.counts.{f}', 'no annotation is read without labels "all"')

        # work, stop, determinism, notes
        w = e['work']
        self.keys(w, ['ranges_visited', 'mask_scans', 'steps']
                  + (['annotation_rows', 'annotation_units', 'memory_bytes'] if labelled else [])
                  + (['extension_edges'] if paths else []), path + '.work')
        self.ok(all(is_int(v) and v >= 0 for v in w.values()), path + '.work')
        self.ok(w['steps'] >= w['ranges_visited'] + w.get('extension_edges', 0), path + '.work',
                'steps >= ranges_visited (+ extension_edges, one step each)')
        stop = e['stop']
        if stop is not None:
            self.keys(stop, SCHEMA['stop'], path + '.stop')
            self.one_of(stop['phase'], STOP_PHASES, path + '.stop.phase')
            self.one_of(stop['reason'], STOP_REASONS, path + '.stop.reason')
            self.ok(total['relation'] != 'exact'
                    or stop['phase'] in ('extraction', 'extension') + ANNOTATION_PHASES,
                    path + '.stop', 'a stop in discovery leaves no exact count')
            self.ok(stop['phase'] not in ANNOTATION_PHASES or labelled, path + '.stop.phase',
                    'an annotation phase without labels "all"')
            self.ok(stop['phase'] != 'extension' or paths, path + '.stop.phase',
                    'an extension without long_search "paths"')
            self.ok((stop['reason'] == 'max_paths') <= (stop['phase'] == 'extension'),
                    path + '.stop.reason', 'max_paths stops the extension only')
            if stop['phase'] == 'extension':
                self.ok(total['relation'] == 'exact'
                        and counts['paths']['relation'] == 'at_least', path + '.stop',
                        'a stop in the extension: the anchors exact, the paths at_least')
            if stop['reason'] == 'max_steps':
                self.ok(w['steps'] <= a['limits']['max_steps'], path + '.work.steps')
        else:
            self.ok(total['relation'] == 'exact', path + '.counts',
                    'a count without a stop is exact')
            if paths:
                self.ok(counts['paths']['relation'] == 'exact'
                        or counts['paths']['extension'] == 'not_admitted', path + '.counts.paths',
                        'paths without a stop: exact, or not admitted')
        self.one_of(e['determinism'], ('full', 'time_limited'), path + '.determinism')
        # SPEC §7.6/§7.9: time_limited iff the clock touched the entry, in its stop or, in
        # partial, only in its cut (stop keeps the first stop: a pattern stopped by max_steps or
        # its threshold whose release met the work time, or a later pattern after a sticky
        # max_steps stop whose empty release met it; review of 2026-10-07, C1-03 and
        # X-DETERMINISM-01). With labels "all", a time stop of the reads after a step stop of the
        # engine states the engine's stop, so only that direction is checked there.
        clocked = (stop is not None and stop['reason'] == 'time') \
            or (e.get('cut') or {}).get('reason') == 'time'
        if clocked:
            self.ok(e['determinism'] == 'time_limited', path + '.determinism',
                    'a time stop or a time cut is not deterministic')
        if e['determinism'] == 'time_limited':
            self.ok(clocked or (labelled and stop is not None), path + '.stop',
                    'time_limited without a time stop or a time cut')
        self.ok(isinstance(e['notes'], list) and all(n in NOTES for n in e['notes']),
                path + '.notes', repr(e['notes']))
        self.ok(e['notes'] == [n for n in NOTES if n in e['notes']], path + '.notes', 'order')
        self.ok(('strand_unknown_canonical' in e['notes']) is (not strand_stated),
                path + '.notes', 'strand_unknown_canonical exactly where no strand is known')
        self.ok(('paths_later_increment' in e['notes']) is (L > k and not paths),
                path + '.notes', 'paths_later_increment exactly for L > k without paths')
        self.ok(('annotation_not_read' in e['notes']) is (named and not labelled),
                path + '.notes', 'annotation_not_read exactly where labels were named, not read')
        if labelled:
            self.ok(('annotation_unbudgeted' in e['notes'])
                    is (e['annotation'] == 'unbudgeted'), path + '.notes')
            self.ok(('record_bounds_unknown' in e['notes']) is (e['placement'] == 'global'),
                    path + '.notes')
            # increment 4: the labels of paths where no coordinate is read
            self.ok(('label_intersection_only' in e['notes'])
                    is (paths and e['placement'] in ('none', 'none_canonical', 'not_requested')),
                    path + '.notes', 'label_intersection_only exactly for paths placed nowhere')
        else:
            self.ok(not {'annotation_unbudgeted', 'record_bounds_unknown',
                         'label_intersection_only'} & set(e['notes']), path + '.notes')
        self.timing(e['timing'], path + '.timing', labelled, paths)

        if mode == 'count':
            self.ok(e['retrieval_complete'] is False, path + '.retrieval_complete',
                    'a count returns no context')
            return
        self.retrieval(e, pk, total, request, a, path, paths)

    def description(self, e, pk, k, path, request):
        strands = request.get('strands', 'both')
        text, L = pk.text, pk.length()
        self.ok(e['pattern'] == text, path + '.pattern', 'the pattern in upper case')
        # L in bases whatever the kind: 3 per residue of a peptide (SPEC §12.2)
        self.ok(e['length'] == L, path + '.length')
        if pk.protein:
            self.ok(e['residues'] == len(text), path + '.residues')
            self.ok(e['genetic_code'] == pk.table, path + '.genetic_code', 'the request\'s')
        self.ok(is_num(e['information_bits']) and abs(e['information_bits'] - pk.bits()) < 1e-9,
                path + '.information_bits', f'{e["information_bits"]}, expected {pk.bits()}')
        if L > k:
            # anchor_information_bits: the bits of P[0, k), whatever the strands (its meaning in
            # contract version 1); min_anchor_information_bits: the least informative anchor
            # window searched, P[0, k) forward, rc(P)[0, k) reverse, whose bits are
            # P[L - k, L)'s (review of 2026-10-07, X-GUARANTEES-01, and the owner's decision);
            # a peptide's windows may cut a codon (exact bits, SPEC §12.2)
            self.ok(is_num(e['anchor_information_bits'])
                    and abs(e['anchor_information_bits'] - pk.bits(0, k)) < 1e-9,
                    path + '.anchor_information_bits', 'the bits of P[0, k)')
            windows = []
            if strands in ('both', 'forward') or pk.palindromic():
                windows.append(pk.bits(0, k))
            if strands in ('both', 'reverse') and not pk.palindromic():
                windows.append(pk.bits(L - k, L))
            self.ok(is_num(e['min_anchor_information_bits'])
                    and abs(e['min_anchor_information_bits'] - min(windows)) < 1e-9,
                    path + '.min_anchor_information_bits', 'the bits of the least anchor window')
        else:
            for f in ('anchor_information_bits', 'min_anchor_information_bits'):
                self.ok(e[f] is None, path + '.' + f)

    @staticmethod
    def exempt(pk, k, request):
        """SPEC §7.8: an exact pattern (every position one base, whatever its kind; a peptide:
        every residue one codon) of L <= k in suffix scope is exempt from the information
        floor."""
        return pk.exact() and pk.length() <= k \
            and request.get('scope', 'any_offset') == 'suffix'

    def retrieval(self, e, pk, total, request, a, path, paths=False):
        k = a['index']['k']
        strand_stated = a['index']['strand_stated']
        mode = a['mode']
        labelled = isinstance(a['output'], dict) and a['output']['labels'] == 'all'
        L = pk.length()
        limits = a['limits']
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
            self.ok(e['cut']['reason'] != 'max_contexts' or L <= k, path + '.cut',
                    'max_contexts cuts contexts')
            self.ok(e['cut']['reason'] != 'max_paths'
                    or (paths and e['returned'] == limits['max_paths']), path + '.cut',
                    'max_paths: the first max_paths paths')
        if e['retrieval_complete']:
            self.ok(total['relation'] == 'exact' and e['stop'] is None
                    and e['withheld'] is None and e['cut'] is None, path,
                    'complete: exact, no stop, nothing withheld or cut')
            if paths:
                c = e['counts']['paths']
                self.ok(c['relation'] == 'exact' and e['returned'] == c['value'], path,
                        'complete: every path returned')
            elif L <= k:
                self.ok(e['returned'] == total['value'], path, 'complete: every context returned')
            else:
                self.ok(e['returned'] == 0 and total['value'] == 0, path,
                        'a long pattern without paths is complete only without anchors')
        elif e['withheld'] is None and mode == 'partial':
            labels_stated = labelled and (
                e['rows_refused'] or e['anchors_truncated'] or e['stop'] is not None
                or e['labels_cut'] is not None or e['occurrences_cut'] is not None)
            self.ok(e['cut'] is not None or labels_stated, path + '.cut',
                    'an incomplete partial list states its cut')
        if mode == 'all_or_count':
            self.ok(e['retrieval_complete'] or e['withheld'] is not None, path,
                    'all_or_count returns all, or withholds')
        if L > k and not paths:
            self.ok(results == [], path + '.results', 'no result without long_search "paths"')
            if not e['retrieval_complete']:
                self.ok(e['withheld'] == {'reason': 'paths_later_increment'}, path + '.withheld')
        if paths:
            self.ok(e['withheld'] != {'reason': 'paths_later_increment'}, path + '.withheld',
                    'paths were asked for')
            c = e['counts']['paths']
            reason = (e['withheld'] or {}).get('reason')
            if reason == 'anchors_above_threshold':
                self.ok(total['relation'] == 'exact' and total['value'] > limits['max_anchors']
                        and c['extension'] == 'not_admitted', path + '.withheld',
                        'the exact anchors above max_anchors, the extension not admitted')
            elif reason == 'count_above_threshold':
                self.ok(c['relation'] == 'exact' and c['value'] > limits['max_paths'],
                        path + '.withheld', 'the exact paths above max_paths')
            elif reason == 'threshold_crossed':
                self.ok(e['stop'] is not None
                        and e['stop']['reason'] in ('max_anchors', 'max_paths'),
                        path + '.withheld', 'a stop_at_threshold stop')
        if mode == 'partial':
            self.ok(e['returned'] <= (limits['max_paths'] if paths else limits['max_contexts']),
                    path + '.returned', 'the cap')

        if paths:
            self.path_results(e, pk, request, a, path)
        else:
            self.context_results(e, pk, a, path)
        if labelled:
            self.labels(e, total, L, a, path, paths, request)

    def path_results(self, e, pk, request, a, path):
        """The paths of long_search "paths" (SPEC §12.1): each with its L bases (sequence =
        instance, an instance of the oriented pattern), its anchor's k bases, offset 0, its node
        path (n = L - k + 1 node ids, the row of each k-mer); never kmer; in the answer order
        (anchor node, orientation, sequence), each once; per strand within the counts."""
        k = a['index']['k']
        strand_stated = a['index']['strand_stated']
        labelled = isinstance(a['output'], dict) and a['output']['labels'] == 'all'
        verified_only = request.get('require_support') == 'record_verified'
        L = pk.length()
        n = L - k + 1
        key = 'strand' if strand_stated else 'orientation'
        allowed = set(e['strands'])
        order = []
        for i, r in enumerate(e['results']):
            rp = f'{path}.results[{i}]'
            self.keys(r, ['sequence', 'anchor_kmer', 'instance', 'offset', key, 'nodes', 'rows']
                      + (['support', 'labels_status', 'labels_total', 'labels'] if labelled
                         else [])
                      + (['labels_excluded_unverified'] if labelled and verified_only else []),
                      rp)
            s = r['sequence']
            self.ok(isinstance(s, str) and len(s) == L and set(s) <= set('ACGT'),
                    rp + '.sequence')
            self.ok(r['instance'] == s, rp + '.instance', 'the sequence')
            self.ok(r['anchor_kmer'] == s[:k], rp + '.anchor_kmer', 'the first k bases')
            self.ok(r['offset'] == 0, rp + '.offset')
            self.ok(r[key] in allowed, rp + '.' + key, f'{r[key]!r} not searched')
            self.ok(pk.instance(s, r[key] in ('-', 'reverse')), rp + '.sequence',
                    f'{s} does not instantiate the oriented pattern')
            nodes, rows = r['nodes'], r['rows']
            self.ok(isinstance(nodes, list) and len(nodes) == n
                    and all(is_int(x) and x >= 1 for x in nodes), rp + '.nodes')
            self.ok(isinstance(rows, list) and len(rows) == n
                    and all(x is None or (is_int(x) and x >= 0) for x in rows), rp + '.rows')
            if a['index']['graph_mode'] == 'basic':
                self.ok(all(x is None or x == y - 1 for x, y in zip(rows, nodes)), rp + '.rows',
                        'node - 1 on a basic graph')
            order.append((nodes[0], ORIENTATION_RANK[r[key]], s))
        self.ok(order == sorted(order) and len(set(order)) == len(order), path + '.results',
                'ordered by (anchor node, orientation, sequence), each path once')
        c = e['counts']['paths']
        by = c['by_strand' if strand_stated else 'by_orientation']
        got = {}
        for r in e['results']:
            part = STRAND_KEYS[r[key]] if strand_stated else r[key]
            got[part] = got.get(part, 0) + 1
        for part, count in by.items():
            m = got.get(part, 0)
            if e['retrieval_complete']:
                self.ok(count['relation'] == 'exact' and m == count['value'],
                        f'{path}.counts.paths.{part}', f'{m} results, the count {count["value"]}')
            elif count['relation'] == 'exact':
                self.ok(m <= count['value'], f'{path}.counts.paths.{part}',
                        f'{m} results above the exact count {count["value"]}')
        self.ok(set(got) <= set(by), path + '.counts.paths', 'a result outside the parts counted')

    def context_results(self, e, pk, a, path):
        k = a['index']['k']
        strand_stated = a['index']['strand_stated']
        labelled = isinstance(a['output'], dict) and a['output']['labels'] == 'all'
        L = pk.length()
        results = e['results']
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
            self.ok(pk.instance(r['instance'], r[key] in ('-', 'reverse')), rp + '.instance',
                    f'{r["instance"]} does not instantiate the oriented pattern')
            self.ok(is_int(r['node']) and r['node'] >= 1, rp + '.node')
            self.ok(r['row'] is None or (is_int(r['row']) and r['row'] >= 0), rp + '.row')
            if a['index']['graph_mode'] == 'basic' and r['row'] is not None:
                self.ok(r['row'] == r['node'] - 1, rp + '.row', 'node - 1 on a basic graph')
            order.append((r['node'], r['offset'], ORIENTATION_RANK[r[key]]))
        self.ok(order == sorted(order) and len(set(order)) == len(order), path + '.results',
                'ordered by (node, offset, orientation), each context once')
        # review of 2026-10-07, X-ORACLE-03: the results and the counts beside them
        contexts = [(r[key], r['kmer'], r['offset']) for r in results]
        self.ok(len(set(contexts)) == len(contexts), path + '.results',
                'a context (strand, k-mer, offset) returned once')
        node_of, kmer_of = {}, {}
        for r in results:
            node_of.setdefault(r['kmer'], set()).add(r['node'])
            kmer_of.setdefault(r['node'], set()).add(r['kmer'])
        self.ok(all(len(v) == 1 for v in node_of.values())
                and all(len(v) == 1 for v in kmer_of.values()), path + '.results',
                'one node per k-mer and one k-mer per node')
        if L <= k and results:
            c = e['counts']['contexts']
            by = c['by_strand' if strand_stated else 'by_orientation']
            per_offset = {}
            per_strand = {}
            for r in results:
                per_offset[str(r['offset'])] = per_offset.get(str(r['offset']), 0) + 1
                s = STRAND_KEYS.get(r[key], r[key]) if strand_stated else r[key]
                per_strand[s] = per_strand.get(s, 0) + 1
            for name, parts, got in (('by_offset', c['by_offset'], per_offset),
                                     ('by_strand' if strand_stated else 'by_orientation', by,
                                      per_strand)):
                for part, count in parts.items():
                    n = got.get(part, 0)
                    if e['retrieval_complete']:
                        # every context returned: the results per part are the part's count
                        self.ok(count['relation'] == 'exact' and n == count['value'],
                                f'{path}.counts.contexts.{name}.{part}',
                                f'{n} results, the count {count["value"]}')
                    elif count['relation'] == 'exact':
                        self.ok(n <= count['value'], f'{path}.counts.contexts.{name}.{part}',
                                f'{n} results above the exact count {count["value"]}')
                self.ok(set(got) <= set(parts), f'{path}.counts.contexts.{name}',
                        'a result outside the parts counted')

    def labels(self, e, total, L, a, path, paths=False, request=None):
        """output.labels "all" (SPEC §14; for paths §12.1): the label fields, their values and the
        rules that tie them to retrieval_complete, withheld and the counts."""
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
        if 'labels_excluded_unverified' in e:
            self.count(e['labels_excluded_unverified'], 'labels',
                       path + '.labels_excluded_unverified')
        if e['withheld'] is not None:
            self.ok(e['by_label'] is None, path + '.by_label', 'null when withheld')
            for f in ('labels', 'occurrences'):
                self.ok(e['counts'][f]['relation'] == 'unknown', f'{path}.counts.{f}',
                        'nothing published when withheld')
            return
        # a path's labels_total is the size of its labels' intersection: unknown (null) when a
        # row of it was truncated (SPEC §12.1)
        unknown_total = ('refused', 'not_read', 'truncated') if paths else ('refused', 'not_read')
        if e['by_label'] is None:
            # partial whose memory account could not hold by_label (SPEC §14.4; outside review
            # GPT-2 recheck): no label of the pattern is built, every context read answers
            # output_budget, the answer is incomplete and states a stop (the output's, or an
            # earlier one: the first stop wins); the counts keep their own relations
            self.ok(mode == 'partial', path + '.by_label', 'null only when withheld or in partial')
            self.ok(not e['retrieval_complete'], path + '.retrieval_complete',
                    'by_label null: not complete')
            self.ok(e['stop'] is not None, path + '.stop', 'by_label null only under a stop')
            self.ok(e['labels_cut'] is None and e['occurrences_cut'] is None, path,
                    'no list is cut when none is built')
            for i, r in enumerate(e['results']):
                rp = f'{path}.results[{i}]'
                self.one_of(r['labels_status'], LABELS_STATUS, rp + '.labels_status')
                self.ok(r['labels_status'] not in ('complete', 'truncated'), rp + '.labels_status',
                        'nothing listed without by_label')
                self.ok(r['labels'] is None, rp + '.labels', 'null without by_label')
                if r['labels_status'] in unknown_total:
                    self.ok(r['labels_total'] is None, rp + '.labels_total')
                elif not paths:
                    self.ok(r['labels_total'] is not None, rp + '.labels_total')
                # a path's excluded count is known only with its labels_total (SPEC §12.1)
                if paths and r.get('labels_excluded_unverified') is not None:
                    self.ok(is_int(r['labels_total']), rp + '.labels_excluded_unverified',
                            'an integer only with the path\'s labels_total')
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
        if paths:
            self.path_labels(e, L, a, path, request, records, placed)
            return
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

    def path_labels(self, e, L, a, path, request, records, placed):
        """The labels of paths (SPEC §12.1, owner decision #14): by_label per label with its
        paths, the paths one record verifies and its occurrences; each path's labels, each with
        its support (record_verified: an occurrence of the whole path in one record, placed;
        label_intersection: on every k-mer of the path, no record holding it whole), the path's
        support theirs; with require_support "record_verified" only the verified ones listed,
        the others counted in labels_excluded_unverified."""
        verified_only = request.get('require_support') == 'record_verified'
        key = 'strand' if a['index']['strand_stated'] else 'orientation'
        keys, verified_labels = [], 0
        for i, b in enumerate(e['by_label']):
            bp = f'{path}.by_label[{i}]'
            self.keys(b, SCHEMA['by_label_paths'], bp)
            self.ok(b['graph'] == a['index']['index_ns'], bp + '.graph')
            self.count(b['paths'], 'paths', bp + '.paths')
            self.count(b['paths_record_verified'], 'paths', bp + '.paths_record_verified')
            self.count(b['occurrences'], 'placed_occurrences', bp + '.occurrences')
            self.ok(b['paths']['value'] >= 1, bp + '.paths', 'a label of a returned path')
            v = b['paths_record_verified']
            if not records:
                self.ok(v['relation'] == 'unknown', bp + '.paths_record_verified',
                        'nothing is verified without a record mapping')
                self.ok(b['occurrences']['relation'] == 'unknown', bp + '.occurrences')
            elif v['relation'] == 'exact' and b['paths']['relation'] == 'exact':
                self.ok(v['value'] <= b['paths']['value'], bp + '.paths_record_verified')
                if verified_only:
                    self.ok(v['value'] == b['paths']['value'], bp + '.paths_record_verified',
                            'require_support: every path listing it verifies it')
            if v['relation'] == 'exact' and v['value'] >= 1:
                verified_labels += 1
            keys.append((-b['paths']['value'], b['column']))
        self.ok(keys == sorted(keys), path + '.by_label', '(paths desc, column asc)')
        support = e['counts']['labels']['by_support']
        if support['record_verified']['relation'] == 'exact' and e['labels_cut'] is None:
            self.ok(support['record_verified']['value'] == verified_labels,
                    path + '.counts.labels.by_support.record_verified',
                    'the labels verified on a returned path')
        if not records:
            self.ok(support['record_verified']['relation'] == 'unknown',
                    path + '.counts.labels.by_support.record_verified')
        rank = {b['column']: i for i, b in enumerate(e['by_label'])}
        for i, r in enumerate(e['results']):
            rp = f'{path}.results[{i}]'
            status = r['labels_status']
            self.one_of(status, LABELS_STATUS, rp + '.labels_status')
            if status == 'complete':
                self.ok(is_int(r['labels_total']), rp + '.labels_total')
            elif status in ('refused', 'not_read', 'truncated'):
                self.ok(r['labels_total'] is None, rp + '.labels_total',
                        'the intersection is known only when every row was read completely')
            if verified_only:
                # the true number or null: an integer only when every row of the path was read
                # completely (as labels_total) and every label carrying it verified or refuted;
                # a path whose verification was not done leaves the entry's count unknown
                excluded = r['labels_excluded_unverified']
                self.ok(excluded is None or (is_int(excluded) and excluded >= 0),
                        rp + '.labels_excluded_unverified')
                if status in ('refused', 'not_read', 'truncated'):
                    self.ok(excluded is None, rp + '.labels_excluded_unverified',
                            'null unless every row of the path was read completely')
                if excluded is not None:
                    self.ok(is_int(r['labels_total']), rp + '.labels_excluded_unverified',
                            'an integer only with the path\'s labels_total')
                elif status == 'complete':
                    self.ok(e['labels_excluded_unverified']['relation'] == 'unknown',
                            rp + '.labels_excluded_unverified',
                            'a path read completely is undecided only when the verification '
                            'was not done, which leaves the entry\'s count unknown')
            if r['labels'] is None:
                self.ok(r['support'] is None, rp + '.support', 'null without labels')
                self.ok(status not in ('complete',), rp + '.labels', 'a complete path lists')
                continue
            self.ok(isinstance(r['labels'], list), rp + '.labels')
            listed = [label['column'] for label in r['labels']]
            self.ok(all(c in rank for c in listed), rp + '.labels', 'labels of by_label')
            self.ok(listed == sorted(listed, key=rank.get), rp + '.labels', 'label order')
            if status == 'complete' and e['labels_cut'] is None:
                if not verified_only:
                    self.ok(len(listed) == r['labels_total'], rp + '.labels',
                            'every label of the intersection')
                elif r['labels_excluded_unverified'] is not None:
                    self.ok(len(listed) + r['labels_excluded_unverified'] == r['labels_total'],
                            rp + '.labels',
                            'every label of the intersection (listed or excluded)')
                else:
                    self.ok(len(listed) <= r['labels_total'], rp + '.labels',
                            'the verified labels of the intersection')
            supports = set()
            for j, label in enumerate(r['labels']):
                lp = f'{rp}.labels[{j}]'
                self.keys(label, ['column', 'support'] + (['occurrences'] if records else [])
                          + (['occurrence_list'] if placed else []), lp)
                self.one_of(label['support'], SUPPORTS, lp + '.support')
                supports.add(label['support'])
                if verified_only:
                    self.ok(label['support'] == 'record_verified', lp + '.support',
                            'require_support: the verified labels only')
                if not records:
                    self.ok(label['support'] == 'label_intersection', lp + '.support',
                            'nothing is verified without a record mapping')
                if records:
                    occ = label['occurrences']
                    self.count(occ, 'placed_occurrences', lp + '.occurrences')
                    if occ['relation'] == 'exact':
                        self.ok((occ['value'] >= 1) is (label['support'] == 'record_verified'),
                                lp + '.support', 'record_verified exactly with an occurrence')
                if not placed or label['occurrence_list'] is None:
                    continue
                self.ok(isinstance(label['occurrence_list'], list), lp + '.occurrence_list')
                if records and label['occurrences']['relation'] == 'exact' \
                        and e['occurrences_cut'] is None:
                    self.ok(len(label['occurrence_list']) == label['occurrences']['value'],
                            lp + '.occurrence_list', 'every occurrence listed')
                for n, o in enumerate(label['occurrence_list']):
                    op = f'{lp}.occurrence_list[{n}]'
                    if records:
                        self.keys(o, SCHEMA['occurrence'], op)
                        m = re.fullmatch(r'(\d+)-(\d+)', o['nt_coords'])
                        self.ok(m and int(m.group(2)) == int(m.group(1)) + L - 1
                                and int(m.group(1)) >= 1
                                and int(m.group(2)) <= o['nt_length'], op + '.nt_coords',
                                'the whole path in the record')
                    else:
                        # a chain of consecutive column coordinates: kmer_coord of its first
                        # k-mer, offset 0 (record bounds unknown: unverified)
                        self.keys(o, SCHEMA['occurrence_global'], op)
                        self.ok(o['offset'] == 0, op + '.offset')
                    self.ok(o['strand'] == r[key], op + '.strand', 'the path\'s')
            expected = None if not supports else supports.pop() if len(supports) == 1 \
                else 'mixed'
            self.ok(r['support'] == expected, rp + '.support',
                    f'{r["support"]!r}, its labels\' {expected!r}')

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
                          ('kinds', KINDS), ('strands', ('both', 'forward', 'reverse')),
                          ('long_search', LONG_SEARCH)):
            self.ok(isinstance(b[f], list) and b[f] and set(b[f]) <= set(values), f'{path}.{f}')
        # increment 4: the paths are opt-in (owner decision #13): an omitted long_search is
        # "anchors", and long_patterns keeps describing that answer
        self.ok(b['default_long_search'] == 'anchors' and 'anchors' in b['long_search'],
                path + '.default_long_search', 'version 1: an omitted long_search is "anchors"')
        self.ok(b['long_patterns'] == 'anchors_counted', path + '.long_patterns')
        # increment 5: the residues, the genetic codes (NCBI ids) and the default
        self.ok(isinstance(b['protein_residues'], list)
                and all(isinstance(r, str) and len(r) == 1 and r in PROTEIN_RESIDUES
                        for r in b['protein_residues']), path + '.protein_residues')
        self.ok(isinstance(b['genetic_codes'], list) and b['genetic_codes']
                and all(c in GENETIC_CODES for c in b['genetic_codes']), path + '.genetic_codes',
                'NCBI translation table ids')
        self.ok(b['default_genetic_code'] in b['genetic_codes'], path + '.default_genetic_code')
        self.ok(isinstance(b['protein_rule'], str) and b['protein_rule'],
                path + '.protein_rule')
        self.ok(('protein' in b['kinds']) is bool(b['protein_residues']), path + '.kinds',
                'protein served exactly with its residues')
        self.ok(b['modes'] == list(MODES), path + '.modes')
        self.ok('none' in b['projections'], path + '.projections', '"none" is served everywhere')
        self.ok(b['default_mode'] in b['modes'], path + '.default_mode')
        self.ok(b['default_projection'] in b['projections'], path + '.default_projection')
        # contract version 1: an omitted output.labels means "none" on every server stating
        # it (review of 2026-10-07, D1-04: a "may become all" would change the meaning of an
        # unchanged request under the same version)
        self.ok(b['default_projection'] == 'none', path + '.default_projection',
                'version 1: an omitted output.labels is "none"')
        # review of 2026-10-07, R2-04: the rule names every cap
        for cap in b['caps']:
            self.ok(cap in b['caps_rule'], path + '.caps_rule', f'{cap} not in the rule')
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
            # SPEC §10.2: a graph the engine does not recognise (review of 2026-10-07, C2-01):
            # not available, and only k is set
            self.one_of(b['unavailable_reason'],
                        ('representation_unsupported', 'primary_unwrapped'),
                        path + '.unavailable_reason')
            for f in ('alphabet', 'strand_stated', 'mask', 'scopes', 'placement', 'support',
                      'annotation'):
                self.ok(b[f] is None, f'{path}.{f}', 'only k is set on an unrecognised graph')
            return
        self.one_of(b['graph_mode'], GRAPH_MODES, path + '.graph_mode')
        self.ok(b['scopes'] == b['scopes_by_graph_mode'][b['graph_mode']], path + '.scopes')
        self.ok(b['strand_stated'] is (b['graph_mode'] == 'basic'), path + '.strand_stated')
        self.one_of(b['mask'], MASKS, path + '.mask')
        self.ok(isinstance(b['alphabet'], str) and b['alphabet'], path + '.alphabet')
        # the alphabet before the mask (pattern.cpp route_support, as the engine orders
        # alphabet_unsupported before mask_required): on the served alphabet the mask is absent
        # exactly when mask_required; on another the alphabet's reason is given, mask or none --
        # a DNA5 graph without a mask is alphabet_untested (review GPT-2 of 2026-10-08,
        # finding 6: this rule held on every alphabet and refused that block)
        if b['alphabet'] == SERVED_ALPHABET:
            self.ok((b['mask'] == 'absent') is (b['unavailable_reason'] == 'mask_required'),
                    path + '.mask', f'on {SERVED_ALPHABET}, mask absent exactly when mask_required')
            self.ok(b['unavailable_reason'] not in ('alphabet_untested', 'alphabet_unsupported'),
                    path + '.unavailable_reason', f'{SERVED_ALPHABET} is served')
        else:
            self.ok(b['available'] is False, path + '.available',
                    f'only {SERVED_ALPHABET} is served')
            self.one_of(b['unavailable_reason'],
                        ('alphabet_untested' if b['alphabet'] == '$ACGTN'
                         else 'alphabet_unsupported',), path + '.unavailable_reason')
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
            check.ok(len(set(p) & set(KINDS)) == 1, f'request.patterns[{i}]', 'one kind')
        check.ok(request.get('long_search', 'anchors') in LONG_SEARCH, 'request.long_search')
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
                                   'slot', 'refusal', 'graph_mode', 'relation', 'note',
                                   'stop_then_cut', 'kind', 'extension', 'support',
                                   'path_support')}
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
                seen['kind'].add(e['kind'])
                if 'error' in e:
                    seen['slot'].add(e['error']['code'])
                    continue
                seen['scope'].add(e['scope'])
                if 'extension' in e['counts'].get('paths', {}):
                    seen['extension'].add(e['counts']['paths']['extension'])
                    for r in e.get('results', []):
                        if 'support' in r:
                            seen['path_support'].add(r['support'])
                            seen['support'].update(x['support'] for x in r['labels'] or [])
                seen['note'].update(e['notes'])
                for f2 in ('withheld', 'cut'):
                    if e.get(f2):
                        seen[f2].add(e[f2]['reason'])
                if e['stop']:
                    seen['stop'].add((e['stop']['phase'], e['stop']['reason']))
                    if e.get('cut') and e['stop']['reason'] != e['cut']['reason']:
                        # an earlier stop kept in stop, a later one shown only by the cut
                        seen['stop_then_cut'].add((e['stop']['reason'], e['cut']['reason'],
                                                   e['determinism']))
                total = e['counts'].get('contexts') or e['counts']['anchors']
                seen['relation'].add(total['relation'])
        self.assertEqual(set(MODES), seen['mode'])
        self.assertEqual(set(SCOPES), seen['scope'])
        # increments 4 and 5: every kind; the paths' extension completed, stopped, not
        # admitted and without anchors; both supports of a label of a path (none listed: null)
        self.assertEqual(set(KINDS), seen['kind'])
        self.assertLessEqual({'completed', 'stopped', 'not_admitted', 'no_anchors'},
                             seen['extension'])
        self.assertEqual(set(SUPPORTS), seen['support'])
        self.assertLessEqual({'record_verified', 'label_intersection', None},
                             seen['path_support'])
        self.assertEqual({'both', 'forward', 'reverse'}, seen['strands'])
        # output_budget from labels_all_output_budget (a GCG repeat whose 1,828 contexts fill
        # the smallest account, max_memory_mb 1); the cut max_memory from labels_all_rows_refused
        # (partial, four patterns sharing that account)
        self.assertEqual(set(WITHHELD), seen['withheld'])
        self.assertEqual({'max_contexts', 'max_steps', 'time', 'max_memory', 'max_paths'},
                         seen['cut'])
        self.assertEqual(set(SLOT_ERRORS), seen['slot'])
        self.assertEqual({'basic', 'primary'}, seen['graph_mode'])
        self.assertEqual(set(NOTES), seen['note'])
        # every relation, bounds and the mask_scan stop included (review of 2026-10-07,
        # X-TESTS-03: no body carried either)
        self.assertEqual(set(RELATIONS), seen['relation'])
        self.assertLessEqual({('discovery', 'max_steps'), ('discovery', 'time'),
                              ('discovery', 'max_contexts'), ('mask_scan', 'max_steps')},
                             seen['stop'])
        # every refusal code but the unproducible ones (review of 2026-10-07, C2-01: the set was
        # pinned to the covered ones, so a code without a fixture went unseen) and the ones no
        # stored fixture holds, listed by name with the index each needs (review GPT-2 of
        # 2026-10-08, finding 6: derived from REFUSALS alone, the expectation could not notice
        # REFUSALS missing two codes; test_the_codes_are_the_sources compares it with the code)
        without_fixture = {'mask_invalid', 'alphabet_untested'}
        self.assertEqual(without_fixture, set(NO_FIXTURE))
        self.assertEqual((set(REFUSALS) | {'initializing'}) - set(UNPRODUCIBLE)
                         - without_fixture, seen['refusal'])
        self.assertEqual(without_fixture | set(UNPRODUCIBLE),
                         set(REFUSALS) - seen['refusal'], 'the refusal codes without a fixture')
        self.assertLessEqual({('label_discovery', 'max_annotation_work')}, seen['stop'])
        # increment 4: the extension's stops
        self.assertLessEqual({('extension', 'max_paths'), ('extension', 'max_steps')},
                             seen['stop'])
        # SPEC §7.6 (the owner's decision of 2026-10-07): stop keeps the first stop; the time
        # that then cut the release shows only as cut time and time_limited, in the pattern
        # stopped by max_steps and in the one after it (review of 2026-10-07, C1-03)
        self.assertIn(('max_steps', 'time', 'time_limited'), seen['stop_then_cut'])

    def test_time_limited_is_the_clocks(self):
        """SPEC §7.6/§7.9, both ways: an entry the clock touched (its stop, or only its cut) is
        time_limited, and a time_limited entry was touched by the clock (review of 2026-10-07,
        C1-03: the validator refused the SPEC's stop max_steps + cut time)."""
        name = 'max_steps_then_time'
        request, answer = self.bodies[name]
        capabilities = self.capabilities_of(self.fixtures[name]['server'])

        class Stub:
            def fail(self, message):
                raise AssertionError(message)

        def check(a):
            Checker(Stub(), name).answer(a, request, capabilities)

        check(answer)
        for i in (0, 1):
            self.assertEqual(('max_steps', 'time', 'time_limited'),
                             (answer['patterns'][i]['stop']['reason'],
                              answer['patterns'][i]['cut']['reason'],
                              answer['patterns'][i]['determinism']))
        # a time cut stated as deterministic
        a = copy.deepcopy(answer)
        for e in a['patterns']:
            e['determinism'] = 'full'
        with self.assertRaisesRegex(AssertionError, 'a time stop or a time cut is not'):
            check(a)
        # time_limited where the clock touched nothing: the first entry's release complete
        # up to max_contexts, cut by its discovery's max_steps
        a = copy.deepcopy(answer)
        a['patterns'][0]['cut'] = {'reason': 'max_steps'}
        with self.assertRaisesRegex(AssertionError, 'time_limited without a time stop'):
            check(a)
        # an entry after a time-limited one stated as deterministic (its own cut dropped)
        a = copy.deepcopy(answer)
        a['patterns'][1]['determinism'] = 'full'
        a['patterns'][1]['cut'] = {'reason': 'max_steps'}
        with self.assertRaisesRegex(AssertionError, 'after a time-limited one'):
            check(a)

    def test_every_unavailable_reason_has_a_capabilities_fixture(self):
        """Every unavailable reason but the unproducible ones and the ones named as without a
        fixture (NO_FIXTURE), on both capabilities routes."""
        seen = {}
        for name, f in self.fixtures.items():
            if f['method'] != 'GET':
                continue
            b = self.bodies[name][1]['pattern']
            if b['available'] is False:
                route = 'probe' if f['path'].startswith('/traverse/') else 'capabilities'
                seen.setdefault(b['unavailable_reason'], set()).add(route)
        without_fixture = {'mask_invalid', 'alphabet_untested'}
        self.assertEqual(without_fixture, set(NO_FIXTURE))
        self.assertEqual(set(UNAVAILABLE) - set(UNPRODUCIBLE) - without_fixture, set(seen))
        self.assertEqual(without_fixture | set(UNPRODUCIBLE), set(UNAVAILABLE) - set(seen),
                         'the unavailable reasons without a fixture')
        for reason, routes in seen.items():
            self.assertEqual({'probe', 'capabilities'}, routes, reason)

    def test_the_codes_are_the_sources(self):
        """REFUSALS and UNAVAILABLE are the codes the server's sources write (review GPT-2 of
        2026-10-08, finding 6: mask_invalid and alphabet_untested, added to pattern.cpp, were
        missing here, and the coverage tests, derived from these lists, could not notice):
        every PatternRefusal's literal code, and every reason support_message words (the
        graph-support reasons a refusal and the capabilities' unavailable_reason carry), with
        the multi-graph server's own unavailable reason."""
        if not os.path.isfile(PATTERN_CPP):
            self.skipTest('the server sources are not in this checkout')
        cli = os.path.dirname(PATTERN_CPP)
        sources = {}
        for name in sorted(os.listdir(cli)):
            if name.endswith('.cpp'):
                with open(os.path.join(cli, name), encoding='utf-8') as f:
                    sources[name] = f.read()
        literal = set()
        for text in sources.values():
            literal |= set(re.findall(r'PatternRefusal\(\s*\d+\s*,\s*"(\w+)"', text))
        support = set(re.findall(r'support\.reason == "(\w+)"', sources['pattern.cpp']))
        self.assertIn('PatternRefusal(400, support.reason, support_message(support))',
                      sources['pattern.cpp'], 'a graph-support refusal is its reason')
        unavailable = set(re.findall(r'p\["unavailable_reason"\] = "(\w+)"',
                                     sources['pattern.cpp']))
        self.assertLessEqual({'mask_required', 'mask_invalid', 'alphabet_untested'}, support,
                             'support_message\'s branches not found in pattern.cpp')
        self.assertEqual(set(REFUSALS), literal | support)
        self.assertEqual(set(UNAVAILABLE), support | unavailable)
        # the slot codes: the engine's PatternError codes (bad_alphabet, stop_unsupported) and
        # its refusals of a parsed pattern (information_below_floor, scope_unsupported)
        engine = os.path.join(REPO, 'src', 'graph', 'alignment', 'pattern_search.cpp')
        with open(engine, encoding='utf-8') as f:
            engine_text = f.read()
        slots = set(re.findall(r'PatternError\(\s*"(\w+)"', engine_text))
        self.assertEqual({'bad_alphabet', 'stop_unsupported'}, slots)
        for code in ('information_below_floor', 'scope_unsupported'):
            self.assertIn(f'"{code}"', engine_text)
        self.assertEqual(set(SLOT_ERRORS), slots | {'information_below_floor',
                                                    'scope_unsupported'})

    def test_by_label_null_in_partial_is_valid(self):
        """SPEC §14.4: in partial, when the memory account cannot hold by_label, the answer has
        by_label null, every context read output_budget and a stop -- a valid v1 answer the
        validator refused (outside review GPT-2 recheck). The body is a real CLI answer on a
        tiny long-label index (data/traverse/pattern_validator/by_label_null_partial); no
        fixture of the mini index can show it (its label names are short)."""
        d = os.path.join(HERE, 'data', 'traverse', 'pattern_validator', 'by_label_null_partial')
        with open(os.path.join(d, 'request.json')) as f:
            request = json.load(f)
        with open(os.path.join(d, 'answer.json')) as f:
            answer = json.load(f)
        e = answer['patterns'][1]
        self.assertEqual((None, {'phase': 'output', 'reason': 'max_memory'}, 'exact'),
                         (e['by_label'], e['stop'], e['counts']['labels']['relation']))

        class Stub:
            def fail(self, message):
                raise AssertionError(message)

        def check(a, req=request):
            Checker(Stub(), 'by_label_null_partial').answer(a, req)

        check(answer)
        # what v1 never answers: a label listed without by_label
        a = copy.deepcopy(answer)
        a['patterns'][1]['results'][0]['labels_status'] = 'complete'
        with self.assertRaisesRegex(AssertionError, 'nothing listed without by_label'):
            check(a)
        # by_label null without a stop (an earlier rule may name it first)
        a = copy.deepcopy(answer)
        a['patterns'][1]['stop'] = None
        with self.assertRaisesRegex(AssertionError, 'by_label null only under a stop|states its cut'):
            check(a)

    def test_the_new_codes_are_valid_answers(self):
        """Hand-made bodies of v1 with mask_invalid and alphabet_untested -- a 400 refusal, and
        the capabilities block of such a server (a DNA5 graph without a mask included) -- are
        accepted (review GPT-2 of 2026-10-08, finding 6: none is stored, see NO_FIXTURE),
        and the rules they rest on still refuse what v1 never answers."""
        class Stub:
            def fail(self, message):
                raise AssertionError(message)

        def refusal(code, status=400):
            check = Checker(Stub(), f'hand-made {code}')
            self.refusal(check, {'status': status, 'headers': {}},
                         {'error': f'pattern: the graph is not served ({code})', 'code': code})

        def block(**fields):
            b = copy.deepcopy(self.capabilities_of('masked'))
            self.assertIsNotNone(b)
            self.assertIs(True, b['available'])
            b.update(fields)
            Checker(Stub(), 'hand-made capabilities').block(b, 'pattern', False)

        for code in ('mask_invalid', 'alphabet_untested'):
            refusal(code)
            with self.assertRaisesRegex(AssertionError, '503 is the deadline'):
                refusal(code, 503)
        # a mask that marks a W = $ edge valid: the graph described, its mask stated
        block(available=False, unavailable_reason='mask_invalid', mask='file')
        # a DNA5 graph, with its mask or without (the alphabet's reason comes first)
        for mask in MASKS:
            block(available=False, unavailable_reason='alphabet_untested', alphabet='$ACGTN',
                  mask=mask)
        # what v1 never answers
        with self.assertRaisesRegex(AssertionError, 'mask absent exactly when mask_required'):
            block(available=False, unavailable_reason='mask_invalid', mask='absent')
        with self.assertRaisesRegex(AssertionError, 'mask absent exactly when mask_required'):
            block(available=False, unavailable_reason='alphabet_untested', mask='absent')
        with self.assertRaisesRegex(AssertionError, 'is served'):
            block(available=False, unavailable_reason='alphabet_untested')
        with self.assertRaisesRegex(AssertionError, 'only \\$ACGT is served'):
            block(alphabet='$ACGTN')
        with self.assertRaisesRegex(AssertionError, 'unavailable_reason'):
            block(available=False, unavailable_reason='mask_required', alphabet='$ACGTN',
                  mask='absent')
        with self.assertRaisesRegex(AssertionError, 'unavailable_reason'):
            block(available=False, unavailable_reason='alphabet_untested', alphabet='$ACGU')

    def test_paths_and_peptides_rules_refuse_what_v1_never_answers(self):
        """Increments 4 and 5 (SPEC §12.1, §12.2): the stored path and peptide answers pass,
        and each rule they rest on refuses an answer that breaks it (so that the rules are not
        vacuous): a path with kmer, a path that does not instantiate the pattern, a label
        stated record_verified without an occurrence, an unverified label listed under
        require_support, a stopped extension stated exact, a peptide's bits, a peptide's
        instance off its codons, a stop refused as bad_alphabet, a path's excluded count stated
        over truncated rows or withheld over complete ones."""
        class Stub:
            def fail(self, message):
                raise AssertionError(message)

        def check(name, mutate=None):
            request, answer = self.bodies[name]
            answer = copy.deepcopy(answer)
            if mutate:
                mutate(answer)
            Checker(Stub(), name).answer(answer, request,
                                         self.capabilities_of(self.fixtures[name]['server']))

        for name in ('paths', 'paths_labels', 'paths_require_support', 'paths_stop_at_max_paths',
                     'peptide', 'peptide_paths', 'peptide_bad_residue', 'paths_global',
                     'paths_primary'):
            check(name)

        def first_path(a):
            return a['patterns'][0]['results'][0]

        cases = [
            ('paths', lambda a: first_path(a).update(kmer=first_path(a)['anchor_kmer']),
             'fields'),
            ('paths', lambda a: first_path(a).update(
                sequence=first_path(a)['sequence'][:-1] + 'ACGT'.replace(
                    first_path(a)['sequence'][-1], '')[0],
                instance=first_path(a)['sequence'][:-1] + 'ACGT'.replace(
                    first_path(a)['sequence'][-1], '')[0]), 'does not instantiate'),
            ('paths_labels',
             lambda a: a['patterns'][1]['results'][0]['labels'][0].update(
                 support='record_verified'), 'record_verified exactly with an occurrence'),
            ('paths_require_support',
             lambda a: a['patterns'][0]['results'][0]['labels'][0].update(
                 support='label_intersection'), 'the verified labels only'),
            ('paths_stop_at_max_paths',
             lambda a: a['patterns'][0]['counts']['paths'].update(relation='exact'),
             'for an extension stopped|an exact total has exact parts'),
            ('paths_stop_at_max_paths',
             lambda a: a['patterns'][0]['counts']['paths'].update(extension='completed'),
             'at_least for an extension completed'),
            ('peptide', lambda a: a['patterns'][0].update(
                information_bits=a['patterns'][0]['information_bits'] + 1), 'information_bits'),
            # the first context's instance: its first codon ATG (M) read as ATA (I in table 1)
            ('peptide', lambda a: a['patterns'][0]['results'][0].update(
                instance='ATA' + a['patterns'][0]['results'][0]['instance'][3:]),
             'instance'),
            ('peptide_bad_residue',
             lambda a: a['patterns'][1]['error'].update(code='bad_alphabet'),
             'bad_alphabet for a pattern inside the alphabet'),
            # a path's labels_excluded_unverified (review of increments 4 and 5, finding 1):
            # an integer for a path whose rows were truncated, and null for a path read and
            # verified completely in an answer whose count is exact
            ('paths_require_support',
             lambda a: first_path(a).update(labels_status='truncated', labels_total=None),
             'null unless every row of the path was read completely'),
            ('paths_require_support',
             lambda a: a['patterns'][1]['results'][0].update(labels_excluded_unverified=None),
             'undecided only when the verification was not done'),
        ]
        for name, mutate, says in cases:
            with self.subTest(fixture=name, says=says):
                with self.assertRaisesRegex(AssertionError, says):
                    check(name, mutate)

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
            # (a multi-graph server and a graph the engine does not recognise state no mask)
            mask = self.bodies[name][1]['pattern'].get('mask') if f['method'] == 'GET' else None
            if mask is not None:
                route = 'probe' if f['path'].startswith('/traverse/') else 'capabilities'
                seen.setdefault(mask, set()).add(route)
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
