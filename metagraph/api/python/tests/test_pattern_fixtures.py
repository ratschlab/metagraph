"""The POST /pattern fixture bodies (data/traverse/pattern/, written by
scripts/traversal/pattern_fixtures.py) against docs/SPEC-pattern-search.md, contract version 1
(with output.labels "all", SPEC §14; the paths of long_search "paths" with their labels'
support, and protein patterns, SPEC §12.1, §12.2; graphs without their dummy-edge mask
answered with bounds and an estimate, their tiny blocks checked and exact, and the stop '*' in
peptides; the work counters, the low-complexity note on completed searches only, the
capabilities' references and their byte budget; and a request's predicate, SPEC §19: the
selection of each pattern of at most k bases, its relations, the projection "predicate_only"
and the predicate block):

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
    by asking a server; a flank of a result's k-mer is checked against the graph by
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
import gzip
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
# output.labels "all" (SPEC §14) adds the annotation phases and reasons; long_search "paths"
# (SPEC §12.1) the phase extension and the reason max_paths; a predicate (SPEC §19.8) the phase
# selection, the reasons max_predicate_work and max_predicate_contexts, three withheld reasons
# and two cut reasons
ANNOTATION_PHASES = ('label_discovery', 'placement', 'output')
STOP_PHASES = ('discovery', 'mask_scan', 'extraction') + ANNOTATION_PHASES + ('extension',
                                                                               'selection')
STOP_REASONS = ('max_steps', 'time', 'max_contexts', 'max_anchors', 'max_annotation_work',
                'max_memory', 'max_paths', 'max_predicate_work', 'max_predicate_contexts')
WITHHELD = ('count_above_threshold', 'threshold_crossed', 'discovery_budget', 'deadline',
            'paths_later_increment', 'annotation_budget', 'anchor_labels_truncated',
            'output_budget', 'anchors_above_threshold', 'predicate_above_threshold',
            'selected_above_threshold', 'predicate_budget')
CUT = ('max_contexts', 'max_steps', 'time', 'max_anchors', 'max_memory', 'max_paths',
       'max_predicate_contexts', 'max_predicate_work',
       # the supported-path search (SPEC §20.5)
       'max_annotation_work')
# in the order an entry lists them (SPEC §7.9): the engine's, then the route's estimate note
# (SPEC §8.10: threshold_upper_bound, no_stop_codon, estimate_sampled_dummy_fraction), then the
# annotation's
NOTES = ('low_complexity_pattern', 'strand_unknown_canonical', 'paths_later_increment',
         'threshold_upper_bound', 'no_stop_codon', 'estimate_sampled_dummy_fraction',
         'annotation_unbudgeted', 'record_bounds_unknown', 'annotation_not_read',
         'label_intersection_only',
         # a predicate (SPEC §19.10)
         'predicate_constant', 'projection_not_read')
# SPEC §7.8: the bases of a pattern its low-complexity diagnostic always reads (one piece of
# sdust, 128 + 63, before any clock reading); a longer one's can be cut by the work time, its
# note then left out
LOW_COMPLEXITY_UNCUT = 191
# (stop_unsupported is a retired slot code: '*' is a residue; no source writes it and no fixture
# holds it)
SLOT_ERRORS = ('bad_alphabet', 'information_below_floor', 'scope_unsupported')
KINDS = ('dna', 'iupac', 'protein')
# what the extension of a pattern longer than k did (counts.paths.extension)
EXTENSIONS = ('no_anchors', 'not_started', 'not_admitted', 'stopped', 'completed')
LONG_SEARCH = ('anchors', 'paths', 'supported_paths')
# long_search "supported_paths" (SPEC §20): the levels a request may ask for, the levels searched
SUPPORTED_PATHS_LEVELS = ('best', 'label_intersection')
# a predicate (SPEC §19): the operators of a predicate, predicate_strands, what a pattern's
# selection did (selection.pass), how it read the predicate's labels (selection.access), what a
# label's presence means (selection.support), the scope of the claim, and the projections
PREDICATE_OPERATORS = ('any', 'all', 'none', 'at_least', 'and', 'or', 'not')
PREDICATE_STRANDS = ('either', 'context')
PASSES = ('completed', 'stopped', 'not_admitted', 'not_started', 'constant')
ACCESS = ('rows', 'columns')
# a selected result's selection_strands (SPEC §19.10)
SELECTION_STRANDS = {'context', 'reverse_complement', 'both', 'either'}
SELECTION_SUPPORTS = ('kmer', 'label_intersection', 'record_verified')
PREDICATE_SCOPES = ('shard_context',)
# the motif-level predicate (SPEC §25): predicate_scope's values, the scope of the claim, what
# decided a motif, why not every context was tested, the evaluation's own stops
PREDICATE_SCOPE_VALUES = ('context', 'motif')
MOTIF_SCOPES = ('shard_motif',)
MOTIF_BASES = ('every_context', 'tested_contexts', 'constant', 'no_instance')
MOTIF_UNTESTED = ('discovery', 'release', 'not_admitted', 'not_started', 'selection')
MOTIF_STOPS = ('time', 'max_memory')
PROJECTIONS = ('none', 'all', 'predicate_only')
# a path result's support: its listed labels' (null when it lists none)
PATH_SUPPORTS = ('record_verified', 'label_intersection', 'mixed')
# peptides: the residues of a protein pattern (any case), the stop '*' among them, and NCBI's
# genetic codes (gc.prt version 4.6, https://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt): each
# table's ncbieaa string, the residue of each codon in TCAG order (TTT TTC TTA TTG TCT ... GGG),
# '*' a stop
PROTEIN_RESIDUES = 'ACDEFGHIKLMNPQRSTVWYXBZJ*'
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
# mask_invalid and alphabet_untested: the route's own reasons not to serve a graph (pattern.cpp
# route_support), a 400 refusal's code and the capabilities' unavailable_reason alike (missing
# here, valid answers carrying them would be rejected)
# mask_required is a code of version 1 that no build serving graphs without their mask answers
# (RETIRED)
REFUSALS = ('invalid_request', 'later_increment', 'resident_only', 'mask_required',
            'representation_unsupported', 'primary_unwrapped', 'alphabet_unsupported', 'deadline',
            'annotation_unbudgeted', 'mask_invalid', 'alphabet_untested',
            # require_support "record_verified" on an index that cannot verify, genetic_code no
            # NCBI table, and a predicate listing more names than caps.max_predicate_labels
            'support_unavailable', 'genetic_code_unknown', 'predicate_too_large')
UNAVAILABLE = ('mask_required', 'representation_unsupported', 'primary_unwrapped',
               'alphabet_unsupported', 'multi_graph_later_increment', 'mask_invalid',
               'alphabet_untested')
# The refusals and unavailable reasons no server of this build gives, so that no fixture holds
# them: primary_unwrapped (server_query and the CLI always wrap a PRIMARY graph in
# CanonicalDBG; only an embedding reaches it) and alphabet_unsupported (the BOSS alphabet is
# the build's: a graph of another alphabet does not load).
UNPRODUCIBLE = ('primary_unwrapped', 'alphabet_unsupported',
                # the block pattern_capabilities_json writes for a multi-graph server that does
                # not serve /pattern: since multi-graph servers serve it (SPEC §24), no route of
                # this build asks for that block
                'multi_graph_later_increment')
# The ones a server of this build gives but no stored fixture holds, each with the index it
# needs: named here rather than left out of the expected sets, so that the coverage tests list
# every code without a fixture. Every other code has a fixture.
NO_FIXTURE = {
    'mask_invalid': 'a graph whose .edgemask marks an edge with W = $ valid: one extended after '
                    'masking by a build older than 85614d30 (metagraph extend on a masked '
                    'graph), or a stale mask',
    'alphabet_untested': 'a DNA5 build of metagraph and a $ACGTN graph (a DNA4 build does not '
                         'load one)',
}
# The codes of version 1 that no build of this kind answers: kept in REFUSALS and UNAVAILABLE,
# since an older build of version 1 answers them and a client keeps handling them, but written
# by no source of this build and held by no fixture.
RETIRED = {
    'mask_required': 'a graph without its dummy-edge mask is served, its counts upper bounds '
                     'with an estimate (capabilities counting upper_bound; the unmasked_* '
                     'fixtures)',
    'resident_only': 'in_ram is accepted, as /search accepts it (the in_ram_single_graph and '
                     'multi_graph_in_ram fixtures)',
}
# the retired codes that were unavailable reasons too (resident_only was a refusal only)
RETIRED_REASONS = {'mask_required'}
# The alphabet the route serves (pattern.cpp alphabet_refusal): on any other the alphabet's
# reason comes first, mask or none
SERVED_ALPHABET = '$ACGT'
MASKS = ('file', 'built_at_load', 'absent')
# how the graph counts (capabilities counting, index.counting): exact with its dummy-edge mask,
# upper_bound without it; the source of the dummy fraction an estimate rests on (SPEC §8.2)
COUNTINGS = ('exact', 'upper_bound')
FRACTION_SOURCES = ('sampled',)
# the units whose bounds carry an estimate on a graph without its mask (the graph's counts)
GRAPH_UNITS = ('graph_contexts', 'anchors', 'paths')
PLACEMENTS = ('record', 'global', 'none', 'none_canonical')
# an entry's placement: the index's, or not_requested (output.occurrences false)
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
                'long_search', 'max_paths', 'require_support', 'genetic_code',
                # a predicate (SPEC §19.2)
                'predicate', 'max_predicate_contexts', 'max_predicate_work', 'predicate_strands',
                # supported paths (SPEC §20.2) and the motif-level predicate (SPEC §25.2)
                'supported_paths_level', 'predicate_scope',
                # the server's (SPEC §24): the graphs of a multi-graph server, the load into RAM
                'graphs', 'in_ram'],
    'request_pattern': ['id', 'dna', 'iupac', 'protein'],
    'request_output': ['labels', 'occurrences', 'paths'],
    'refusal': ['error', 'code'],
    'answer': ['pattern_contract_version', 'mode', 'output', 'index', 'limits', 'timing',
               'patterns', 'predicate'],
    'index': ['index_ns', 'index_fp', 'release', 'k', 'graph_mode', 'alphabet',
              'strand_stated', 'counting', 'dummy_fraction'],
    # SPEC §8.2: the dummy fraction of a graph without its mask
    'dummy_fraction': ['value', 'interval', 'samples', 'source'],
    'limits': ['max_contexts', 'max_anchors', 'max_steps', 'time_budget_ms',
               'finalize_reserve_ms', 'min_information_bits', 'max_patterns',
               'stop_at_threshold', 'max_labels_per_anchor', 'max_annotation_work',
               'max_memory_mb', 'max_labels', 'max_occurrences_per_label',
               'allow_unbudgeted_annotation', 'long_search', 'max_paths', 'require_support',
               'max_predicate_contexts', 'max_predicate_work', 'max_predicate_labels',
               'predicate_strands',
               # supported paths (SPEC §20.2), the motif-level predicate (SPEC §25.2)
               'supported_paths_level', 'predicate_scope', 'clamped'],
    'clamped': ['field', 'requested', 'effective'],
    'timing': ['elapsed_ms', 'label_discovery_ms', 'placement_ms', 'extension_ms',
               # the labels of paths (SPEC §12.1)
               'label_intersection_ms', 'verification_ms',
               # a predicate (SPEC §19.10)
               'selection_ms',
               # the supported-path search (SPEC §20.6)
               'support_ms',
               # in_ram (SPEC §24)
               'load_ms'],
    'entry': ['id', 'kind', 'pattern', 'length', 'residues', 'genetic_code', 'information_bits',
              'anchor_information_bits', 'min_anchor_information_bits', 'error', 'mode', 'scope',
              'strands', 'palindromic', 'counts', 'work', 'stop',
              'retrieval_complete', 'withheld', 'returned', 'cut', 'results', 'absence_scope',
              'determinism', 'notes', 'timing', 'placement', 'annotation', 'by_label',
              'rows_refused', 'anchors_truncated', 'labels_cut', 'occurrences_cut',
              'labels_excluded_unverified',
              # a predicate (SPEC §19.10)
              'selection', 'absence_filter',
              # the motif-level predicate (SPEC §25.3)
              'motif'],
    'counts': ['contexts', 'anchors', 'paths', 'labels', 'occurrences', 'tested', 'selected',
               # long_search "supported_paths" (SPEC §20.4)
               'supported_paths'],
    'count': ['value', 'relation', 'unit', 'lower', 'upper', 'estimate'],
    'contexts_count': ['suffix', 'by_offset', 'by_strand', 'by_orientation'],
    'anchors_count': ['by_strand', 'by_orientation'],
    # long_search "paths" (SPEC §12.1)
    'paths_count': ['by_strand', 'by_orientation', 'candidates_examined', 'extension'],
    # long_search "supported_paths" (SPEC §20.4)
    'supported_paths_count': ['by_strand', 'by_orientation', 'level', 'search',
                              'candidates_examined', 'branches_pruned',
                              'branches_pruned_by_predicate'],
    'labels_count': ['by_support'],
    'by_support': ['record_verified', 'label_intersection'],
    'work': ['ranges_visited', 'mask_scans', 'steps', 'annotation_rows', 'annotation_units',
             'memory_bytes', 'extension_edges',
             # counters beside the steps (SPEC §8.7)
             'extension_anchors', 'extension_branches', 'annotation_rows_distinct',
             'verification_steps',
             # a predicate (SPEC §19.10)
             'predicate_rows', 'predicate_units', 'predicate_lookups',
             # the supported-path search (SPEC §20.6, §20.9), the motif (SPEC §25.3)
             'anchor_rows', 'row_cache_hits', 'row_cache_evictions', 'mirror_rows',
             'motif_units'],
    'stop': ['phase', 'reason'],
    'reason': ['reason'],
    'result': ['kmer', 'instance', 'offset', 'strand', 'orientation', 'node', 'row', 'support',
               'labels_status', 'labels_total', 'labels', 'selection_labels',
               'selection_strands'],
    # a path (a result of a pattern longer than k under long_search "paths")
    'path_result': ['sequence', 'anchor_kmer', 'instance', 'offset', 'strand', 'orientation',
                    'nodes', 'rows', 'support', 'labels_status', 'labels_total', 'labels',
                    'labels_excluded_unverified',
                    # a selected supported path (SPEC §20.9)
                    'selection_labels', 'selection_strands'],
    'error': ['code', 'message'],
    'capabilities': ['pattern_contract_version', 'available', 'unavailable_reason', 'modes',
                     'default_mode', 'projections', 'default_projection',
                     'projections_later_increment', 'kinds', 'kinds_later_increment',
                     'default_scope', 'scopes_by_graph_mode', 'scopes', 'long_patterns',
                     'strands', 'default_strands', 'graph_cleaned', 'records_shorter_than_k',
                     'resident_only', 'caps', 'default_time_budget_ms', 'finalize_reserve_ms',
                     'caps_rule', 'graph_mode', 'k', 'alphabet', 'strand_stated', 'mask',
                     'placement', 'support', 'annotation', 'default_occurrences',
                     # paths and peptides
                     'long_search', 'default_long_search', 'protein_residues', 'genetic_codes',
                     'default_genetic_code', 'protein_rule',
                     # a graph without its mask
                     'counting', 'dummy_fraction',
                     # the delivery rates and the prose references (SPEC §10.2)
                     'delivery_mbps',
                     # a predicate (SPEC §19.12)
                     'predicate',
                     # the max_anchors and max_labels of a request that names none, at most
                     # their caps (SPEC §4.5)
                     'default_max_anchors', 'default_max_labels',
                     # in_ram accepted (SPEC §24)
                     'in_ram'],
    'delivery_mbps': ['build', 'compress'],
    # a predicate (SPEC §19): the predicate's operators (each a one-member object), the
    # answer's predicate block, an entry's selection, and the capabilities' predicate object
    'predicate': list(PREDICATE_OPERATORS),
    'predicate_block': ['normal_form', 'names', 'known', 'unknown_labels', 'vacuous', 'scope',
                        'strands', 'motif_scope'],
    'selection': ['pass', 'support', 'access'],
    'capabilities_predicate': ['operators', 'strands', 'scopes', 'access'],
    # the motif-level predicate (SPEC §25.3): an entry's motif block and a label found
    'motif': ['selected', 'decided_by', 'untested', 'labels', 'labels_present', 'labels_absent',
              'stop'],
    'motif_label': ['column', 'contexts', 'strands'],
    # a multi-graph server (SPEC §24): its answer, the tags of each pair's answer, the block of
    # one pair (GET /pattern/capabilities?graph=), and GET /capabilities' graph_summary
    'answer_multi': ['pattern_contract_version', 'graphs', 'answered', 'refused', 'answers',
                     'timing'],
    'answer_pair': ['graph', 'graph_path', 'annotation_path', 'index_fp', 'outcome', 'refusal'],
    # a refused entry's `refusal`: the body the pair alone would have been answered, with its
    # HTTP status (no code for an unexpected failure)
    'pair_refusal': ['http_status', 'error', 'code'],
    'capabilities_pair': ['graph', 'graph_path'],
    'graph_summary': ['columns_disjoint', 'shared_columns', 'pairs'],
    'graph_summary_pair': ['graph', 'graph_path', 'annotation_path', 'index_ns', 'index_fp', 'k',
                           'graph_mode', 'available', 'unavailable_reason', 'mask', 'counting',
                           'traversal'],
    'graph_summary_traversal': ['regime', 'num_labels', 'has_coordinates',
                                'has_coord_to_header', 'supports_trace'],
    # the block of GET /traverse/capabilities (SPEC §23): the gate fields and the route of the
    # full block
    'pattern_gate': None,
    # output.labels "all" (SPEC §14)
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
# The fields of the full block that GET /traverse/capabilities keeps for a client that reads it
# alone (SPEC §23, the gate block): every key the search service gates on or parses, and
# counting -- its interface inventory of 2026-10-08, §3.10: the version, available and
# unavailable_reason gate the host (capabilities.py:637-645, 752-781); modes, projections,
# kinds, protein_residues, genetic_codes and strands gate request values (validate.py:667,
# 676-684, 338-343, 147-175, 520-531, 688); scopes, scopes_by_graph_mode, graph_mode and k gate
# scopes and lengths (capabilities.py:694, 728-736, 809-814); long_patterns and long_search the
# long patterns (capabilities.py:738-750); finalize_reserve_ms is the budget floor
# (validate.py:561-575); support, placement and annotation gate record_verified and labels
# (validate.py:485-515, 685-687); default_long_search, default_genetic_code,
# default_occurrences, default_time_budget_ms and mask are parsed (capabilities.py:670-714).
# The server's list is pattern.cpp's kPatternGateKeys and kPatternGateCaps
GATE_KEYS = ['pattern_contract_version', 'available', 'unavailable_reason',
             'modes', 'projections', 'kinds', 'protein_residues', 'genetic_codes', 'strands',
             'scopes', 'scopes_by_graph_mode', 'graph_mode', 'k',
             'long_patterns', 'long_search', 'default_long_search',
             'finalize_reserve_ms', 'default_time_budget_ms', 'default_genetic_code',
             'default_occurrences',
             'support', 'placement', 'annotation', 'mask', 'counting',
             'caps']
# the caps of the gate block: the ceilings and the chunk size (max_patterns) the service gates
# on (validate.py:423-460, 701-717), and the two it parses (min_information_bits, max_anchors)
GATE_CAPS = ['max_patterns', 'max_contexts', 'max_anchors', 'max_paths', 'max_steps',
             'time_budget_ms', 'min_information_bits', 'max_memory_mb', 'max_labels_per_anchor']
SCHEMA['pattern_gate'] = GATE_KEYS + ['details']
# the route of the full block, which the /traverse block names (SPEC §23)
DETAILS = 'GET /pattern/capabilities'
# what caps_rule refers to: the section that classifies every cap (SPEC §4.5)
CAPS_RULE = 'SPEC-pattern-search.md section 4.5'
# the full block's graph fields, null while the index loads (SPEC §10.2)
GRAPH_FIELDS = ('graph_mode', 'k', 'alphabet', 'strand_stated', 'mask', 'counting',
                'dummy_fraction', 'scopes', 'placement', 'support', 'annotation')
ENTRY_DESCRIPTION = ['id', 'kind', 'pattern', 'length', 'information_bits',
                     'anchor_information_bits', 'min_anchor_information_bits']
ENTRY_ANSWERED = ENTRY_DESCRIPTION + ['mode', 'scope', 'strands', 'palindromic', 'counts',
                                      'work', 'stop', 'retrieval_complete', 'absence_scope',
                                      'determinism', 'notes', 'timing']
ENTRY_RETRIEVAL = ['withheld', 'returned', 'cut', 'results']
# output.labels "all" in a retrieval mode
ENTRY_LABELS = ['placement', 'annotation', 'by_label', 'rows_refused', 'anchors_truncated',
                'labels_cut', 'occurrences_cut']
LIMITS_LABELS = ['max_labels_per_anchor', 'max_annotation_work', 'max_memory_mb', 'max_labels',
                 'max_occurrences_per_label']
# (max_checked_entries: the server's limit of the unchecked candidates a pattern on a graph
# without its mask has tested; no request field)
CAPS = ['max_contexts', 'max_anchors', 'max_steps', 'time_budget_ms', 'min_information_bits',
        'max_patterns'] + LIMITS_LABELS + ['max_paths', 'max_checked_entries',
                                           # a predicate (SPEC §19.12)
                                           'max_predicate_contexts', 'max_predicate_work',
                                           'max_predicate_labels',
                                           # the row cache of a supported-path search (SPEC
                                           # §20.3): the server's ceiling, no request field
                                           'row_cache_mb']
# the limits a predicate answer echoes (SPEC §19.10)
LIMITS_PREDICATE = ['max_predicate_contexts', 'max_predicate_work', 'max_predicate_labels',
                    'predicate_strands']


def load(*path):
    with open(os.path.join(DATA, *path), encoding='utf-8') as f:
        return json.load(f)


def served_bytes(value):
    """The length of |value| as the server writes it (jsoncpp, compact JSON, every byte outside
    ASCII escaped as \\uXXXX), from above: each float counted as 24 characters (jsoncpp writes
    up to 17 significant digits, json.dumps the shortest that reads back), the rest exactly."""
    if isinstance(value, float):
        return 24
    # the brackets, the commas between the members, a colon per key
    if isinstance(value, dict):
        return 2 + max(len(value) - 1, 0) + sum(served_bytes(k) + 1 + served_bytes(v)
                                                for k, v in value.items())
    if isinstance(value, list):
        return 2 + max(len(value) - 1, 0) + sum(served_bytes(v) for v in value)
    return len(json.dumps(value, ensure_ascii=True))


# The route's sources (src/cli/*.cpp) that the route does not call, with the refusal codes they
# write (a source moves out of this list once the route calls it, its codes into REFUSALS)
NOT_SERVED_SOURCES = {}

# The ceiling of the capabilities document the service's MCP tool returns in one piece
# (api/python/metagraph/traverse/mcp_tools.py, CAPABILITIES_MAX_BYTES; read from its source by
# test_capabilities_documents_keep_a_kibibyte, no import), and the room every fixture server's
# document keeps under it
MCP_TOOLS = os.path.join(REPO, 'api', 'python', 'metagraph', 'traverse', 'mcp_tools.py')
CAPABILITIES_MAX_BYTES = 32 * 1024
CAPABILITIES_ROOM = 1024


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


def predicate_names(p):
    """The names a predicate (or a normal form) lists, in the order of their first appearance."""
    out = []

    def walk(q):
        if isinstance(q, bool):
            return
        (op, v), = q.items()
        if op in ('any', 'all', 'none'):
            names = v
        elif op == 'at_least':
            names = v['labels']
        elif op == 'not':
            walk(v)
            return
        else:
            for x in v:
                walk(x)
            return
        for n in names:
            if n not in out:
                out.append(n)
    walk(p)
    return out


def predicate_holds(p, labels):
    """A predicate (or a normal form, true and false included) on a set of column names: the
    validator's own evaluator (SPEC §19.3), never the server's."""
    if isinstance(p, bool):
        return p
    (op, v), = p.items()
    if op in ('any', 'all', 'none'):
        hits = sum(n in labels for n in v)
        return hits > 0 if op == 'any' else hits == len(v) if op == 'all' else hits == 0
    if op == 'at_least':
        return sum(n in labels for n in v['labels']) >= v['n']
    if op == 'and':
        return all(predicate_holds(q, labels) for q in v)
    if op == 'or':
        return any(predicate_holds(q, labels) for q in v)
    return not predicate_holds(v, labels)


def predicate_kleene(p, present):
    """A normal form's strong-Kleene value (True, False, or None: undecided) when the labels of
    |present| are present for sure and every other label of it may or may not be (SPEC §25.4):
    definite exactly when every completion gives the same value (the validator's own
    evaluator, by its definition: every completion of the undecided labels tried)."""
    names = [n for n in predicate_names(p) if n not in present]
    values = set()
    for bits in range(1 << len(names)):
        extra = {n for i, n in enumerate(names) if bits >> i & 1}
        values.add(predicate_holds(p, set(present) | extra))
        if len(values) == 2:
            return None
    return values.pop()


def predicate_fold(p, unknown):
    """SPEC §19.4: the normal form of a request's predicate once the names of |unknown| are folded
    away (the validator's own folding)."""
    (op, v), = p.items()
    if op in ('any', 'all', 'none', 'at_least'):
        listed = v['labels'] if op == 'at_least' else v
        known = [n for n in listed if n not in unknown]
        if op == 'any' and not known:
            return False
        if op == 'all' and len(known) < len(listed):
            return False
        if op == 'none' and not known:
            return True
        if op == 'at_least':
            if len(known) < v['n']:
                return False
            return {'at_least': {'n': v['n'], 'labels': known}}
        return {op: known}
    if op == 'not':
        q = predicate_fold(v, unknown)
        return (not q) if isinstance(q, bool) else {'not': q}
    absorbing = op == 'or'
    kept = []
    for q in v:
        f = predicate_fold(q, unknown)
        if isinstance(f, bool):
            if f == absorbing:
                return absorbing
            continue
        kept.append(f)
    if not kept:
        return not absorbing
    if len(kept) == 1:
        return kept[0]
    return {op: kept}


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
    B: D or N, Z: E or Q, J: I or L; '*': the stop codons -- none in the tables without an
    unconditional stop codon, 27, 28 and 31, whose context stops are their residue's; a stop
    codon never matches any other residue)."""
    if residue == '*':
        return {c for c, a in zip(CODONS_TCAG, GENETIC_CODES[table]) if a == '*'}
    admits = {'X': lambda a: True, 'B': lambda a: a in 'DN', 'Z': lambda a: a in 'EQ',
              'J': lambda a: a in 'IL'}.get(residue, lambda a: a == residue)
    return {c for c, a in zip(CODONS_TCAG, GENETIC_CODES[table]) if a != '*' and admits(a)}


def round_half_up(x):
    """std::round for x >= 0 (half away from zero), not Python's round (half to even)."""
    r = math.floor(x)
    return r + 1 if x - r >= 0.5 else r


def estimate_of(count, fraction):
    """SPEC §7.4: a bounds count's estimate on a graph without its mask, round(upper x f) kept
    inside [lower, upper]."""
    return max(count['lower'], min(count['upper'], round_half_up(count['upper'] * fraction)))


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
                # a residue admitting no codon ('*' in a table without a stop codon) counts as
                # one exact codon, so that the bits stay finite (SPEC §12.2)
                spelled = len({c[lo:hi] for c in codons}) or 1
                out += 2 * (hi - lo) - math.log2(spelled)
        return out

    def has_instances(self):
        """False iff no string instantiates the pattern: a peptide holding '*' read in a table
        without an unconditional stop codon (27, 28, 31)."""
        return not self.protein or all(self.codons)

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

    def answered_unsearched(self, k):
        """SPEC §12.2: a pattern without instances of at most k bases is answered exact 0
        without a search, before the information floor and the budget (note no_stop_codon)."""
        return not self.has_instances() and self.length() <= k


class Checker:
    """Validates one answer body; every failure names the fixture and the JSON path."""

    def __init__(self, test, name):
        self.t = test
        self.name = name
        # the dummy fraction of a graph served without its mask (the answer's
        # index.dummy_fraction.value), None with the mask; set by answer()
        self.fraction = None
        # whether a count of the entry being checked states an estimate
        self.estimated = False
        # on a graph without its mask whose capabilities are known, the server's
        # max_checked_entries; and the answer's graph mode. Set by answer()
        self.checked_limit = None
        self.graph_mode = None
        # the normal form of the answer's predicate (None without one, or when its binding
        # stopped), and the capabilities of the fixture's server; set by answer()
        self.normal_form = None
        self.capabilities = None

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

    def count(self, c, unit, path, extra=(), graph=True):
        keys = ['value', 'relation', 'unit'] + list(extra)
        # SPEC §7.4: on a graph without its mask every bounds count of the graph (contexts,
        # anchors, paths) carries its estimate, and no other count does (not the selection's
        # tested and selected, |graph| false: SPEC §19.7)
        estimate = isinstance(c, dict) and c.get('relation') == 'bounds' \
            and self.fraction is not None and unit in GRAPH_UNITS and graph
        if isinstance(c, dict) and c.get('relation') == 'bounds':
            keys += ['lower', 'upper']
        if estimate:
            keys += ['estimate']
        self.keys(c, keys, path)
        if estimate:
            self.estimated = True
            self.ok(is_int(c['estimate']) and c['estimate'] == estimate_of(c, self.fraction),
                    path + '.estimate', f'{c["estimate"]!r}, expected round(upper x f) in '
                    f'[lower, upper] = {estimate_of(c, self.fraction)}')
        self.one_of(c['relation'], RELATIONS, path + '.relation')
        self.ok(c['unit'] == unit, path + '.unit', f'{c["unit"]!r}, expected {unit!r}')
        if c['relation'] == 'unknown':
            self.ok(c['value'] is None, path + '.value', 'unknown counts are null')
        else:
            self.ok(is_int(c['value']) and c['value'] >= 0, path + '.value', repr(c['value']))
        if c['relation'] == 'bounds':
            self.ok(is_int(c['lower']) and is_int(c['upper']) and c['lower'] <= c['upper']
                    and c['value'] == c['lower'], path, 'bounds: value == lower <= upper')

    def dummy_fraction(self, f, path):
        """SPEC §8.2: f, the fraction of real k-mers among the entries a pattern can count, with
        its 95% interval (Wilson's, holding f), the entries drawn and its source."""
        self.keys(f, SCHEMA['dummy_fraction'], path)
        self.ok(is_num(f['value']) and 0 <= f['value'] <= 1, path + '.value')
        self.ok(isinstance(f['interval'], list) and len(f['interval']) == 2
                and all(is_num(x) for x in f['interval'])
                and 0 <= f['interval'][0] <= f['value'] <= f['interval'][1] <= 1,
                path + '.interval', 'an interval in [0, 1] holding the value')
        self.ok(is_int(f['samples']) and f['samples'] > 0, path + '.samples')
        self.one_of(f['source'], FRACTION_SOURCES, path + '.source')
        if is_int(f['samples']) and f['samples'] > 0 and is_num(f['value']):
            # value = real / samples: a whole number of the entries drawn
            real = f['value'] * f['samples']
            self.ok(abs(real - round(real)) < 1e-6, path + '.value', 'real / samples')

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

    # ------------------------------------------------------------- predicates (5b)

    def predicate_form(self, p, path, depth=0):
        """A predicate written in the request's syntax (SPEC §19.3): one operator per object,
        lists of distinct non-empty strings, at_least's n within its list, nesting bounded."""
        self.ok(isinstance(p, dict) and len(p) == 1, path, 'one operator per predicate')
        if not isinstance(p, dict) or len(p) != 1:
            return
        (op, v), = p.items()
        self.one_of(op, PREDICATE_OPERATORS, path)
        self.ok(depth <= 64, path, 'nested at most 64 deep')
        if op in ('any', 'all', 'none') or op == 'at_least':
            names = v['labels'] if op == 'at_least' else v
            if op == 'at_least':
                self.keys(v, ['n', 'labels'], f'{path}.at_least')
                self.ok(is_int(v['n']) and 1 <= v['n'] <= len(names), f'{path}.at_least.n')
            self.ok(isinstance(names, list) and names
                    and all(isinstance(n, str) and n for n in names)
                    and len(set(names)) == len(names), f'{path}.{op}',
                    'a list of distinct non-empty names')
        elif op == 'not':
            self.predicate_form(v, f'{path}.not', depth + 1)
        else:
            self.ok(isinstance(v, list) and v, f'{path}.{op}', 'a list of predicates')
            for i, q in enumerate(v):
                self.predicate_form(q, f'{path}.{op}[{i}]', depth + 1)

    def predicate_block(self, b, request, index, path):
        """The answer's predicate block (SPEC §19.10): the request's predicate as bound to this
        index -- its normal form (the validator's own folding of the request's predicate with the
        unknown names, SPEC §19.4), the names, the known ones, the unknown ones in the order of
        their first appearance, vacuous (the normal form on a context carrying none of its
        labels), the scope, the strands evaluated (either on CANONICAL and PRIMARY graphs). A
        binding the account or the time stopped states the names only (the rest null)."""
        motif = request.get('predicate_scope') == 'motif'
        self.keys(b, [f for f in SCHEMA['predicate_block'] if motif or f != 'motif_scope'], path)
        if motif:
            # SPEC §25.3: the scope of each pattern's motif claim
            self.one_of(b['motif_scope'], MOTIF_SCOPES, path + '.motif_scope')
        asked = request['predicate']
        names = predicate_names(asked)
        self.ok(b['names'] == len(names), path + '.names', 'the distinct names of the request')
        self.one_of(b['scope'], PREDICATE_SCOPES, path + '.scope')
        self.one_of(b['strands'], PREDICATE_STRANDS, path + '.strands')
        expected = request.get('predicate_strands', 'either') \
            if index['graph_mode'] == 'basic' else 'either'
        self.ok(b['strands'] == expected, path + '.strands',
                f'{b["strands"]!r}, evaluated {expected!r}')
        if b['normal_form'] is None:
            for f in ('known', 'unknown_labels', 'vacuous'):
                self.ok(b[f] is None, f'{path}.{f}', 'null with the normal form')
            return
        unknown = b['unknown_labels']
        self.ok(isinstance(unknown, list) and all(isinstance(n, str) for n in unknown)
                and unknown == [n for n in names if n in unknown], path + '.unknown_labels',
                'names of the request, in the order of their first appearance')
        self.ok(b['known'] == len(names) - len(unknown), path + '.known', 'names - unknown')
        nf = b['normal_form']
        if not isinstance(nf, bool):
            self.predicate_form(nf, path + '.normal_form')
        self.ok(nf == predicate_fold(asked, set(unknown)), path + '.normal_form',
                f'{nf!r}, the request folded: {predicate_fold(asked, set(unknown))!r}')
        self.ok(b['vacuous'] is predicate_holds(nf, set()), path + '.vacuous',
                'the normal form on a context carrying none of its labels')
        self.normal_form = nf

    def selection(self, e, total, L, k, request, a, path):
        """An entry's selection (SPEC §19.7): selection {pass, support, access}, counts.tested and
        counts.selected with the relations of the pass, selected <= tested <= the raw count, the
        work of a pass that read nothing."""
        s = e['selection']
        self.keys(s, SCHEMA['selection'], path + '.selection')
        self.one_of(s['pass'], PASSES, path + '.selection.pass')
        self.one_of(s['access'], ACCESS, path + '.selection.access')
        self.one_of(s['support'], SELECTION_SUPPORTS, path + '.selection.support')
        self.ok(e['absence_filter'] == 'predicate', path + '.absence_filter')
        unit = 'graph_contexts' if L <= k else 'paths'
        T, S = e['counts']['tested'], e['counts']['selected']
        self.count(T, unit, path + '.counts.tested', graph=False)
        self.count(S, unit, path + '.counts.selected', graph=False)
        # a pattern longer than k: its supported paths selected (SPEC §20.9), the counts'
        # relations those of §19.7 with the supported paths as the raw count; under
        # long_search "anchors" never started (§19.5)
        supported = L > k and request.get('long_search') == 'supported_paths'
        if supported:
            self.ok(s['support'] == e['counts']['supported_paths']['level'],
                    path + '.selection.support', 'the level the supported paths were searched at')
            self.ok(e['work']['predicate_rows'] == 0, path + '.work.predicate_rows',
                    'supported paths: the selection reads no row of its own')
        elif L > k:
            self.ok(s['pass'] == 'not_started', path + '.selection.pass',
                    'a pattern longer than k: not_started')
        else:
            self.ok(s['support'] == 'kmer', path + '.selection.support', 'L <= k: kmer')
        nf = self.normal_form
        bound = isinstance(a.get('predicate'), dict) and a['predicate']['normal_form'] is not None
        if L <= k or supported:
            self.ok((s['pass'] == 'constant') is (bound and isinstance(nf, bool)),
                    path + '.selection.pass', 'constant exactly for a constant normal form')
        if not bound:
            self.ok(s['pass'] == 'not_started', path + '.selection.pass',
                    'a predicate not bound: no pass')
        p = s['pass']
        R = e['counts']['supported_paths'] if supported else total
        r_upper = R['upper'] if R['relation'] == 'bounds' else R['value']
        if p == 'completed':
            # (supported paths a monotone predicate pruned are not counted: at_least, every
            # one completed decided)
            self.ok(T['relation'] == 'exact'
                    and R['relation'] in (('exact', 'at_least') if supported else ('exact',))
                    and T['value'] == R['value'], path + '.counts.tested',
                    'completed: every raw context tested')
            self.ok(S['relation'] == 'exact', path + '.counts.selected', 'completed: exact')
        elif p == 'stopped':
            self.ok(T['relation'] == 'exact', path + '.counts.tested', 'the decisions made')
            if R['relation'] in ('exact', 'bounds'):
                self.ok(S['relation'] == 'bounds'
                        and S['upper'] == S['lower'] + r_upper - T['value'],
                        path + '.counts.selected', 'stopped: bounds [S, S + R_upper - T]')
            else:
                self.ok(S['relation'] == 'at_least', path + '.counts.selected',
                        'stopped on a raw count without an upper bound: at_least')
        elif p == 'not_admitted':
            self.ok(T['relation'] == 'unknown' and S['relation'] == 'unknown', path + '.counts',
                    'not_admitted: unknown')
        elif p == 'not_started':
            self.ok(T['relation'] == 'unknown', path + '.counts.tested', 'not_started: unknown')
            no_raw = (R['relation'], R['value']) == ('exact', 0)
            self.ok(S == ({'value': 0, 'relation': 'exact', 'unit': unit} if no_raw
                          else {'value': None, 'relation': 'unknown', 'unit': unit}),
                    path + '.counts.selected', 'not_started: unknown, exact 0 without a context')
        else:
            # constant: tested the raw count; false selects nothing, true the raw count
            raw = {f: v for f, v in R.items() if f in ('value', 'relation', 'unit', 'lower',
                                                        'upper')}
            self.ok(T == raw, path + '.counts.tested', 'constant: the raw count')
            self.ok(S == (raw if nf is True else {'value': 0, 'relation': 'exact', 'unit': unit}),
                    path + '.counts.selected', 'constant: 0, or the raw count')
        if T['relation'] != 'unknown' and S['relation'] != 'unknown':
            self.ok(S['value'] <= T['value'], path + '.counts', 'selected <= tested')
        if T['relation'] != 'unknown' and R['relation'] in ('exact', 'bounds'):
            self.ok(T['value'] <= r_upper, path + '.counts', 'tested <= the raw count')
        w = e['work']
        if p in ('constant', 'not_admitted', 'not_started'):
            self.ok(w['predicate_rows'] == 0, path + '.work.predicate_rows', f'{p}: no row read')
        if p in ('constant', 'not_admitted'):
            self.ok(w['predicate_units'] == 0 and w['predicate_lookups'] == 0, path + '.work',
                    f'{p}: no work')
        self.ok(w['predicate_lookups'] == 0 or (a['index']['graph_mode'] == 'basic'
                                                 and a['predicate']['strands'] == 'either'),
                path + '.work.predicate_lookups', 'reverse-complement lookups: either on BASIC')

    def supported_paths_count(self, counts, anchors, e, strand_stated, request, a, predicate,
                              path):
        """counts.supported_paths and counts.paths of a pattern longer than k under long_search
        "supported_paths" (SPEC §20.4): the supported paths with their split, the level searched
        (the one asked for, else the index's best), what the search did, the branches it entered
        and pruned (with a predicate those it pruned, which leave the supported paths at_least);
        counts.paths a plain count, never below the supported paths."""
        c = counts['supported_paths']
        sp = path + '.supported_paths'
        key = 'by_strand' if strand_stated else 'by_orientation'
        self.count(c, 'paths', sp, extra=[key, 'level', 'search', 'candidates_examined',
                                          'branches_pruned']
                   + (['branches_pruned_by_predicate'] if predicate else []))
        self.by_orientation(c, e, strand_stated, 'paths', sp)
        self.one_of(c['level'], SUPPORTS, sp + '.level')
        if request.get('supported_paths_level') == 'label_intersection':
            self.ok(c['level'] == 'label_intersection', sp + '.level', 'the level asked for')
        elif self.capabilities is not None:
            self.ok(c['level'] == self.capabilities['support'], sp + '.level',
                    'the index\'s best support')
        self.one_of(c['search'], EXTENSIONS, sp + '.search')
        for f in ('candidates_examined', 'branches_pruned') \
                + (('branches_pruned_by_predicate',) if predicate else ()):
            self.ok(is_int(c[f]) and c[f] >= 0, f'{sp}.{f}')
        pruned = c.get('branches_pruned_by_predicate', 0)
        relation = {'completed': 'at_least' if pruned else 'exact', 'no_anchors': 'exact',
                    'stopped': 'at_least', 'not_started': 'unknown',
                    'not_admitted': 'unknown'}.get(c['search'])
        self.ok(c['relation'] == relation, sp + '.relation',
                f'{c["relation"]} for a search {c["search"]}'
                + (' that a predicate pruned' if pruned else ''))
        no_anchor = (anchors['relation'], anchors['value']) == ('exact', 0)
        self.ok((c['search'] == 'no_anchors') is no_anchor, sp + '.search',
                'no_anchors exactly when the anchors are exact 0')
        if no_anchor:
            self.ok(c['value'] == 0 and c['candidates_examined'] == 0, sp,
                    'no anchor, no path: exact 0')
        if c['search'] == 'not_admitted':
            self.ok(self.above(anchors, a['limits']['max_anchors']), sp + '.search',
                    'not admitted: the exact anchors (or their upper bound) above max_anchors')
        if c['relation'] in ('exact', 'at_least'):
            self.ok(c['value'] <= c['candidates_examined'] or c['value'] == 0,
                    sp + '.candidates_examined', 'at least one branch per path')
        if predicate and pruned:
            # only a monotone normal form under "context" prunes (SPEC §20.9)
            self.ok(a['predicate']['strands'] == 'context'
                    or a['index']['graph_mode'] != 'basic', sp + '.branches_pruned_by_predicate',
                    'with "either" the predicate does not prune')
        # counts.paths: a plain count of the complete walks (SPEC §20.4)
        p = counts['paths']
        self.ok(set(p) <= {'value', 'relation', 'unit', 'lower', 'upper', 'estimate'},
                path + '.paths', 'a plain count beside the supported paths')
        if c['search'] in ('not_started', 'not_admitted'):
            self.ok(p['relation'] == 'unknown', path + '.paths', 'the extension did not run')
        elif c['search'] == 'stopped':
            self.ok(p['relation'] == 'at_least', path + '.paths', 'a stop in the search')
        elif no_anchor:
            self.ok((p['relation'], p['value']) == ('exact', 0), path + '.paths')
        elif c['branches_pruned'] == 0 and not pruned:
            self.ok(p['relation'] == 'exact', path + '.paths',
                    'a completed search that pruned nothing counts every walk')
        if p['relation'] in ('exact', 'at_least') and c['relation'] in ('exact', 'at_least'):
            self.ok(c['value'] <= p['value'], path + '.paths',
                    'every supported path is a complete walk')
        if pruned:
            self.ok(p['relation'] == 'at_least', path + '.paths',
                    'walks a predicate pruned are not counted')

    def motif(self, e, total, L, k, a, path):
        """An entry's motif block (SPEC §25.3, §25.4): the normal form on the union of its
        contexts' labels -- decided by every context (labels_present every label found, the
        value the normal form on them, labels_absent the others), by the labels found on the
        tested ones (a presence claim only: Kleene-definite on them), by a constant normal form,
        by the pattern having no instance (the empty, complete union: no_instance, never
        every_context), or undecided with the cause; labels present in label order, each with
        its contexts and strands; a pattern longer than k undecided, but one without anchors
        (no_instance)."""
        m = e['motif']
        mp = path + '.motif'
        self.keys(m, SCHEMA['motif'], mp)
        nf = self.normal_form
        names = predicate_names(nf) if isinstance(nf, dict) else []
        self.ok(m['labels'] == len(names), mp + '.labels', 'the labels of the normal form')
        basis = m['decided_by']
        self.ok(basis is None or basis in MOTIF_BASES, mp + '.decided_by')
        self.ok(m['untested'] is None or m['untested'] in MOTIF_UNTESTED, mp + '.untested')
        self.ok(m['stop'] is None or m['stop'] in MOTIF_STOPS, mp + '.stop')
        self.ok((m['selected'] is None) is (basis is None), mp + '.selected',
                'a value exactly when decided')
        self.ok(m['selected'] in (None, True, False), mp + '.selected')
        present = m['labels_present']
        if present is not None:
            self.ok(isinstance(present, list), mp + '.labels_present')
            keys = []
            for i, l in enumerate(present):
                lp = f'{mp}.labels_present[{i}]'
                self.keys(l, SCHEMA['motif_label'], lp)
                self.ok(l['column'] in names, lp + '.column', 'a label of the normal form')
                self.count(l['contexts'], 'graph_contexts', lp + '.contexts', graph=False)
                self.ok(l['contexts']['relation']
                        == ('exact' if basis in ('every_context', 'no_instance') else 'at_least'),
                        lp + '.contexts', 'exact with every context tested, else at_least')
                self.ok(is_int(l['contexts']['value']) and l['contexts']['value'] >= 1,
                        lp + '.contexts', 'found on a context')
                self.one_of(l['strands'], SELECTION_STRANDS, lp + '.strands')
                self.ok((l['strands'] == 'either') is (a['index']['graph_mode'] != 'basic'),
                        lp + '.strands', '"either" exactly where one row serves both')
                if a['index']['graph_mode'] == 'basic' \
                        and a['predicate']['strands'] == 'context':
                    self.ok(l['strands'] == 'context', lp + '.strands',
                            'predicate.strands "context": the own row only')
                keys.append((-l['contexts']['value'], l['column']))
                if e['counts'].get('tested', {}).get('relation') == 'exact':
                    self.ok(l['contexts']['value'] <= e['counts']['tested']['value'],
                            lp + '.contexts', 'at most the contexts tested')
            self.ok(keys == sorted(keys) and len(set(keys)) == len(keys),
                    mp + '.labels_present', 'label order (contexts desc, column asc)')
        found = {l['column'] for l in present or []}
        no_context = (total['relation'], total['value']) == ('exact', 0)
        if basis == 'every_context':
            self.ok(m['untested'] is None and present is not None, mp,
                    'every context tested: the union listed')
            self.ok(m['selected'] is predicate_holds(nf, found), mp + '.selected',
                    'the normal form on the union')
            self.ok(m['labels_absent'] == m['labels'] - len(found), mp + '.labels_absent',
                    'the labels of the normal form not found')
            self.ok(L <= k and e['selection']['pass'] == 'completed' and not no_context,
                    mp + '.decided_by',
                    'every context tested: the pass completed on at least one context')
        elif basis == 'no_instance':
            self.ok(no_context, mp + '.decided_by', 'no instance: the raw count exact 0')
            self.ok(m['untested'] is None and m['stop'] is None and present == [], mp,
                    'no instance: nothing to test, nothing found')
            self.ok(m['selected'] is predicate_holds(nf, set()), mp + '.selected',
                    'the normal form on the empty set')
            self.ok(m['labels_absent'] == m['labels'], mp + '.labels_absent',
                    'every label of the normal form')
        elif basis == 'tested_contexts':
            self.ok(m['untested'] is not None and present is not None, mp,
                    'decided by the labels found: some context untested')
            self.ok(m['labels_absent'] is None, mp + '.labels_absent',
                    'no absence claim from an incomplete pass')
            self.ok(m['selected'] is predicate_kleene(nf, found), mp + '.selected',
                    'definite on the labels found whatever the others')
        elif basis == 'constant':
            self.ok(isinstance(nf, bool) and m['selected'] is nf, mp + '.selected',
                    'a constant normal form')
            self.ok(m['untested'] is None and present == [] and m['labels_absent'] == 0, mp)
        else:
            self.ok(m['labels_absent'] is None, mp + '.labels_absent', 'undecided: no claim')
            self.ok(m['untested'] is not None or m['stop'] is not None or nf is None, mp,
                    'undecided: a context untested, or a stop')
            if isinstance(nf, dict) and present is not None and m['untested'] is not None:
                self.ok(predicate_kleene(nf, found) is None or m['stop'] is not None,
                        mp + '.selected', 'undecided although the labels found decide it')
        if L > k and not no_context and basis != 'constant':
            self.ok(m['untested'] == 'not_started' and basis is None, mp,
                    'a pattern longer than k: not asked of its walks')
        if no_context and isinstance(nf, dict) and m['stop'] is None:
            self.ok(basis == 'no_instance', mp + '.decided_by',
                    'a pattern without an instance says so (constant first)')
        if m['stop'] is not None:
            self.ok(e['stop'] is not None, path + '.stop', 'the motif\'s stop on the entry')
        if m['stop'] == 'time':
            self.ok(present is None and basis is None, mp, 'a time stop: undecided, unlisted')
        self.ok(e['work']['motif_units'] <= e['work']['predicate_units'],
                path + '.work.motif_units', 'part of predicate_units')

    def rows_refused(self, e, path):
        """The rows the account refused (SPEC §14.4): each with its phase (with a predicate:
        the selection's, phase selection, first)."""
        phases = [x.get('phase') for x in e['rows_refused']]
        for i, x in enumerate(e['rows_refused']):
            xp = f'{path}.rows_refused[{i}]'
            self.keys(x, SCHEMA['row_refused'], xp)
            # (extension: the supported-path search's own reads, SPEC §20.5)
            self.one_of(x['phase'], ('label_discovery', 'placement', 'selection', 'extension'),
                        xp + '.phase')
            self.ok(x['reason'] == 'max_memory', xp + '.reason')
        self.ok(phases == sorted(phases, key=lambda p: p != 'selection'), path + '.rows_refused',
                'the selection\'s first')

    def selection_labels(self, e, path, a=None, level=None):
        """The results' selection_labels: the predicate's labels in the set each context was
        evaluated on, names of the normal form, satisfying it, in label order (contexts desc,
        column asc over the results); beside them selection_strands: per label the orientation
        whose row carries it -- "either" on CANONICAL and PRIMARY graphs, "context" for every
        label under predicate.strands "context", and otherwise "context" or "both" exactly when
        the label is on the context's own row (its labels, when the projection read them) and
        "reverse_complement" when it is not."""
        names = set(predicate_names(self.normal_form)) \
            if isinstance(self.normal_form, dict) else set()
        basic = a is None or a['index']['graph_mode'] == 'basic'
        either = a is None or a['predicate']['strands'] == 'either'
        lists = []
        for i, r in enumerate(e['results']):
            rp = f'{path}.results[{i}].selection_labels'
            got = r.get('selection_labels')
            self.ok(isinstance(got, list), f'{path}.results[{i}]',
                    f'fields: selection_labels missing')
            if not isinstance(got, list):
                continue
            self.ok(all(isinstance(n, str) for n in got) and len(set(got)) == len(got), rp,
                    'a list of distinct names')
            self.ok(set(got) <= names, rp, 'names of the normal form')
            self.ok(predicate_holds(self.normal_form, set(got)), rp,
                    'a selected context\'s labels satisfy the normal form')
            lists.append(got)
            sp = f'{path}.results[{i}].selection_strands'
            strands = r.get('selection_strands')
            self.ok(isinstance(strands, list) and len(strands) == len(got),
                    f'{path}.results[{i}]', 'fields: selection_strands, one per selection label')
            if not isinstance(strands, list) or len(strands) != len(got):
                continue
            self.ok(all(x in SELECTION_STRANDS for x in strands), sp,
                    f'one of {sorted(SELECTION_STRANDS)}')
            if not basic:
                self.ok(all(x == 'either' for x in strands), sp,
                        'one row for both orientations: "either"')
                continue
            self.ok('either' not in strands, sp, 'a BASIC graph states the row')
            if not either:
                self.ok(all(x == 'context' for x in strands), sp,
                        'predicate.strands "context": the own row only')
            own = r.get('labels')
            if isinstance(own, list):
                # (a supported path's own support at the record level: its verified labels)
                own = {label['column'] for label in own
                       if level != 'record_verified' or label['support'] == 'record_verified'}
                for n, x in zip(got, strands):
                    self.ok((n in own) == (x in ('context', 'both')), sp,
                            f'{n!r} {x}: on the own row exactly when "context" or "both"')
        count = {}
        for got in lists:
            for n in got:
                count[n] = count.get(n, 0) + 1
        for i, got in enumerate(lists):
            self.ok(got == sorted(got, key=lambda n: (-count[n], n)),
                    f'{path}.results[{i}].selection_labels', 'label order')

    # ------------------------------------------------------------- answers

    def answer(self, a, request, capabilities=None):
        self.keys(a, [f for f in SCHEMA['answer'] if f != 'predicate' or 'predicate' in request],
                  'answer')
        self.ok(a['pattern_contract_version'] == 1, 'pattern_contract_version')
        mode = request.get('mode', 'all_or_count')
        self.ok(a['mode'] == mode, 'mode', f'{a["mode"]!r}, the request asks {mode!r}')
        asked = request.get('output', {})
        # a predicate request, and the projections that read labels ("all", and
        # "predicate_only" with a predicate)
        predicate = 'predicate' in request
        labelled = mode != 'count' and asked.get('labels') in ('all', 'predicate_only')
        if mode == 'count':
            self.ok(a['output'] is None, 'output', 'null in mode count')
        elif labelled:
            self.ok(a['output'] == {'labels': asked['labels'],
                                    'occurrences': asked.get('occurrences', True)},
                    'output', repr(a['output']))
        else:
            self.ok(a['output'] == {'labels': 'none'}, 'output', repr(a['output']))
        self.ok(('predicate' in a) is predicate, 'predicate',
                'the predicate block exactly with a predicate request')
        self.normal_form = None
        self.capabilities = capabilities
        if predicate:
            self.predicate_block(a['predicate'], request, a['index'], 'predicate')

        index = a['index']
        # SPEC §8.2: counting and dummy_fraction only in the answers on a graph without its
        # mask; an answer on a masked graph has neither
        unmasked = 'counting' in index
        self.keys(index, [f for f in SCHEMA['index']
                          if unmasked or f not in ('counting', 'dummy_fraction')], 'index')
        for f in ('index_ns', 'index_fp'):
            self.ok(index[f] is None or isinstance(index[f], str), 'index.' + f)
        self.ok(isinstance(index['release'], str), 'index.release')
        self.ok(is_int(index['k']) and index['k'] >= 2, 'index.k')
        self.one_of(index['graph_mode'], GRAPH_MODES, 'index.graph_mode')
        self.one_of(index['alphabet'], ('$ACGT', '$ACGTN'), 'index.alphabet')
        self.ok(index['strand_stated'] is (index['graph_mode'] == 'basic'),
                'index.strand_stated', 'true exactly on basic graphs')
        self.fraction = None
        if unmasked:
            self.ok(index['counting'] == 'upper_bound', 'index.counting',
                    'stated only on a graph without its mask')
            self.dummy_fraction(index['dummy_fraction'], 'index.dummy_fraction')
            self.fraction = index['dummy_fraction']['value']
        self.graph_mode = index['graph_mode']
        self.checked_limit = capabilities['caps']['max_checked_entries'] \
            if capabilities is not None and unmasked else None
        if capabilities is not None:
            # the graph's counting and its dummy fraction, sampled once: the capabilities'
            self.ok(capabilities['counting'] == ('upper_bound' if unmasked else 'exact'),
                    'index', 'the capabilities\' counting')
            self.ok(index.get('dummy_fraction') == capabilities['dummy_fraction'],
                    'index.dummy_fraction', 'the capabilities\' dummy_fraction')

        limits = a['limits']
        # the annotation limits only in the answers that read annotation; the path limits only
        # in the answers to long_search "paths" (SPEC §12.1), require_support among them only
        # when labels are read
        paths = request.get('long_search', 'anchors') == 'paths'
        # long_search "supported_paths" (SPEC §20.2): its search reads rows in every mode (the
        # annotation limits echoed), with the level asked for
        supported = request.get('long_search') == 'supported_paths'
        # with a predicate the account is used in every mode (the annotation limits echoed),
        # and the selection's limits are echoed (predicate_scope when the request named it)
        annotated = labelled or predicate or supported
        self.keys(limits, [f for f in SCHEMA['limits']
                           if (annotated
                               or f not in LIMITS_LABELS + ['allow_unbudgeted_annotation'])
                           and (paths or supported or f not in ('long_search', 'max_paths'))
                           and ((paths or supported) and labelled or f != 'require_support')
                           and (supported or f != 'supported_paths_level')
                           and (predicate or f not in LIMITS_PREDICATE)
                           and ('predicate_scope' in request or f != 'predicate_scope')],
                  'limits')
        if 'predicate_scope' in request:
            self.ok(limits['predicate_scope'] == request['predicate_scope'],
                    'limits.predicate_scope', 'as requested')
        if supported:
            self.ok(limits['long_search'] == 'supported_paths', 'limits.long_search')
            self.ok(limits['supported_paths_level']
                    == request.get('supported_paths_level', 'best'),
                    'limits.supported_paths_level', 'as requested, default applied')
            self.ok(is_int(limits['max_paths']) and limits['max_paths'] >= 0, 'limits.max_paths')
            if labelled:
                self.ok(limits['require_support']
                        == request.get('require_support', 'label_intersection'),
                        'limits.require_support', 'as requested')
        if predicate:
            for f in ('max_predicate_contexts', 'max_predicate_work', 'max_predicate_labels'):
                self.ok(is_int(limits[f]) and limits[f] >= 0, 'limits.' + f)
            self.ok(limits['max_predicate_work'] >= 1, 'limits.max_predicate_work')
            self.ok(limits['predicate_strands'] == request.get('predicate_strands', 'either'),
                    'limits.predicate_strands', 'as requested')
            if capabilities is not None:
                self.ok(limits['max_predicate_labels'] == capabilities['caps']['max_predicate_labels'],
                        'limits.max_predicate_labels', 'the server\'s')
        if annotated:
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
            + ['max_paths', 'max_predicate_contexts', 'max_predicate_work']
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

        # in_ram (SPEC §24): the time before the work began, waiting for and loading the index
        self.timing(a['timing'], 'timing', load='in_ram' in request)
        self.ok(len(a['patterns']) == len(request['patterns']), 'patterns',
                'one entry per request pattern')
        for i, (entry, asked) in enumerate(zip(a['patterns'], request['patterns'])):
            self.entry(entry, asked, request, a, f'patterns[{i}]')
        self.budget(a, request, limits)

    def multi_answer(self, a, request, pairs):
        """A multi-graph server's answer (SPEC §24.1): the names answered, one entry per pair
        (|pairs| per name on the fixture server) each tagged with its pair and its index_fp and
        stating its outcome -- answered (the pair's answer) or refused (`refusal`: the body
        the pair alone would have been answered, with its HTTP status, and nothing of an
        answer) -- the counts of the outcomes, and the request's time"""
        self.keys(a, SCHEMA['answer_multi'], 'answer')
        self.ok(a['pattern_contract_version'] == 1, 'pattern_contract_version')
        asked = request.get('graphs')
        self.ok(a['graphs'] == sorted(set(asked)) if asked is not None
                else a['graphs'] == sorted(a['graphs']), 'graphs',
                'the names, deduplicated, in byte order')
        self.ok(len(a['answers']) == pairs * len(a['graphs']), 'answers', 'one per pair')
        self.keys(a['timing'], ['elapsed_ms'], 'timing')
        names = [x['graph'] for x in a['answers']]
        self.ok(names == sorted(names) and set(names) == set(a['graphs']), 'answers',
                'in the order of the names')
        tags = [f for f in SCHEMA['answer_pair'] if f != 'refusal']
        outcomes = [x.get('outcome') for x in a['answers']]
        self.ok(all(o in ('answered', 'refused') for o in outcomes), 'answers',
                'each entry\'s outcome, answered or refused')
        self.ok(a['answered'] == outcomes.count('answered')
                and a['refused'] == outcomes.count('refused'), 'answered',
                'answered and refused count the entries of each outcome')
        for i, x in enumerate(a['answers']):
            path = f'answers[{i}]'
            self.ok(set(tags) <= set(x), path, 'tagged with its pair')
            self.ok(isinstance(x['graph_path'], str) and isinstance(x['annotation_path'], str),
                    path + '.graph_path')
            if x['outcome'] == 'refused':
                self.keys(x, SCHEMA['answer_pair'], path)
                self.ok(x['index_fp'] is None or isinstance(x['index_fp'], str),
                        path + '.index_fp')
                self.pair_refusal(x['refusal'], path + '.refusal')
                continue
            self.ok('refusal' not in x, path, 'an answered entry carries no refusal')
            self.ok(x['index_fp'] == x['index']['index_fp'], path + '.index_fp',
                    'the pair\'s index_fp')

    def pair_refusal(self, r, path):
        """The refusal of one pair of a multi-graph request (SPEC §24.1): the body the pair
        alone would have been answered (a refusal's {error, code}; a failure's {error}) with
        its HTTP status: 503 is the deadline, every refusal with a code else 400, a failure
        without a code 400 or 500."""
        self.ok(isinstance(r, dict) and set(r) <= set(SCHEMA['pair_refusal'])
                and {'http_status', 'error'} <= set(r), path, 'http_status, error[, code]')
        self.ok(isinstance(r['error'], str) and r['error'], path + '.error')
        if 'code' in r:
            self.one_of(r['code'], REFUSALS, path + '.code')
            self.ok((r['http_status'] == 503) is (r['code'] == 'deadline'), path + '.http_status',
                    '503 is the deadline, every other refusal 400')
            self.ok(r['http_status'] in (400, 503), path + '.http_status')
        else:
            self.ok(r['http_status'] in (400, 500), path + '.http_status',
                    'a failure without a code: 400 or 500')

    def budget(self, a, request, limits):
        """SPEC §7.6, one budget per request: the steps of all patterns sum to at most max_steps,
        and a stop by max_steps or time is sticky -- every later answered pattern stops with the
        same reason in discovery, its counts unknown, its (engine) work zero. Except a pattern
        answered without a search (a peptide without instances of at most k bases, note
        no_stop_codon, SPEC §12.2): exact 0, no stop, work zero, determinism full, wherever it
        sits, like an error slot."""
        answered = [e for e in a['patterns'] if 'error' not in e]
        self.ok(sum(e['work']['steps'] for e in answered) <= limits['max_steps'], 'patterns',
                'the steps of all patterns sum to at most max_steps')
        sticky = None
        timed = False
        k = a['index']['k']
        for i, e in enumerate(a['patterns']):
            if 'error' in e:
                continue
            path = f'patterns[{i}]'
            if Kind(request['patterns'][i], request).answered_unsearched(k):
                # (memory_bytes, with labels "all", is the request's account so far)
                engine = ('ranges_visited', 'mask_scans', 'steps', 'extension_edges',
                          'extension_anchors', 'extension_branches')
                self.ok(e['stop'] is None and e['determinism'] == 'full'
                        and all(e['work'][f] == 0 for f in engine if f in e['work']),
                        path, 'answered without a search: no stop, determinism full, no work')
                continue
            # SPEC §7.9: the entry a time stop touched and every answered entry after it
            if timed:
                self.ok(e['determinism'] == 'time_limited', path + '.determinism',
                        'an entry after a time-limited one is time-limited')
            timed = timed or e['determinism'] == 'time_limited'
            if sticky is not None:
                self.ok(e['stop'] == {'phase': 'discovery', 'reason': sticky}, path + '.stop',
                        f'the {sticky} stop of an earlier pattern, in discovery')
                engine = ('ranges_visited', 'mask_scans', 'steps', 'extension_edges',
                          'extension_anchors', 'extension_branches')
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

    def timing(self, t, path, labelled=False, extension=False, selection=False, load=False,
               supported=False):
        # (the supported-path search, SPEC §20.6: its extension and its reads, no label read
        # of their own)
        if supported:
            labelled, extension = False, False
        self.keys(t, ['elapsed_ms'] + (['load_ms'] if load else [])
                  + (['label_discovery_ms', 'placement_ms'] if labelled else [])
                  + (['extension_ms'] if extension or supported else [])
                  + (['support_ms'] if supported else [])
                  + (['label_intersection_ms', 'verification_ms'] if labelled and extension
                     else [])
                  + (['selection_ms'] if selection else []), path)
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
            # (without the mask, the anchors' upper bound above max_anchors)
            self.ok(self.above(anchors, limits['max_anchors']), path + '.extension',
                    'not admitted: the exact anchors (or their upper bound) above max_anchors')
        if c['relation'] in ('exact', 'at_least'):
            # each path entered its last branch of its own
            self.ok(c['value'] <= c['candidates_examined'], path + '.candidates_examined',
                    'at least one branch per path')

    def above(self, c, threshold):
        """A threshold decided against a count with no stop: the exact count above it, or on a
        graph without its mask (conservative) a bounds count whose upper bound is above it."""
        if c['relation'] == 'exact':
            return c['value'] > threshold
        return self.fraction is not None and c['relation'] == 'bounds' \
            and c['upper'] > threshold

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
        labelled = isinstance(a['output'], dict) and a['output']['labels'] in ('all',
                                                                               'predicate_only')
        # a predicate request (its selection in every answered entry)
        predicate = 'predicate' in request
        # the request named the labels (or an annotation field) and this answer reads none
        named = (request.get('output', {}).get('labels') == 'all'
                 or any(f in request for f in LIMITS_LABELS + ['allow_unbudgeted_annotation',
                                                               'require_support']))
        # a field only a projection reads (SPEC §19.10, projection_not_read); for a pattern
        # searched for its supported paths, whose search reads under max_annotation_work and
        # never reads max_labels_per_anchor, fewer (SPEC §20.6)
        projection_named = (request.get('output', {}).get('labels') in ('all', 'predicate_only')
                            or any(f in request for f in ('max_labels_per_anchor',
                                                          'max_annotation_work', 'max_labels',
                                                          'max_occurrences_per_label',
                                                          'require_support')))
        supported_projection_named = (
            request.get('output', {}).get('labels') in ('all', 'predicate_only')
            or any(f in request for f in ('max_labels', 'max_occurrences_per_label',
                                          'require_support')))
        # predicate_scope "motif" (SPEC §25): every answered entry states its motif
        motif = predicate and request.get('predicate_scope') == 'motif'
        self.ok(set(e) <= set(SCHEMA['entry']), path, f'fields outside the SPEC: '
                f'{sorted(set(e) - set(SCHEMA["entry"]))}')
        self.ok(e['id'] == asked.get('id'), path + '.id', 'echoed')
        pk = Kind(asked, request)
        self.ok(e['kind'] == pk.kind, path + '.kind')
        text = pk.text
        # a peptide's entry states its residues and its genetic code
        description = ENTRY_DESCRIPTION + (['residues', 'genetic_code'] if pk.protein else [])

        self.estimated = False
        if 'error' in e:
            err = e['error']
            self.keys(err, SCHEMA['error'], path + '.error')
            self.one_of(err['code'], SLOT_ERRORS, path + '.error.code')
            self.ok(isinstance(err['message'], str) and err['message'], path + '.error.message')
            if err['code'] == 'bad_alphabet':
                # (the stop '*' is a residue, SPEC §12.2)
                self.keys(e, ['id', 'kind', 'error'], path)
                self.ok(not pk.parsed(), path, 'bad_alphabet for a pattern inside the alphabet')
                return
            # SPEC §12.2: a pattern without instances of at most k bases is answered, exact 0,
            # before the information floor (no slot error)
            self.ok(not pk.answered_unsearched(k), path,
                    'a slot error for a pattern without instances of at most k bases')
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
        # long_search "paths": a pattern longer than k answered by its paths; "supported_paths":
        # by its supported paths (SPEC §20)
        paths = request.get('long_search', 'anchors') == 'paths' and L > k
        supported = request.get('long_search') == 'supported_paths' and L > k
        verified_only = request.get('require_support') == 'record_verified'
        expected = description + ENTRY_ANSWERED[len(ENTRY_DESCRIPTION):] \
            + (ENTRY_RETRIEVAL if mode != 'count' else []) + (ENTRY_LABELS if labelled else []) \
            + (['labels_excluded_unverified'] if labelled and (paths or supported)
               and verified_only else []) \
            + (['selection', 'absence_filter'] if predicate else []) \
            + (['rows_refused'] if (predicate or supported) and not labelled else []) \
            + (['motif'] if motif else [])
        self.keys(e, expected, path)
        self.description(e, pk, k, path, request)
        # SPEC §7.8: a pattern answered is exempt or at or above the floor (for L > k: every
        # searched anchor window)
        bits = e['min_anchor_information_bits'] if L > k else e['information_bits']
        self.ok(self.exempt(pk, k, request) or bits >= a['limits']['min_information_bits']
                or pk.answered_unsearched(k), path, 'answered below the floor')
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
        selection_counts = ['tested', 'selected'] if predicate else []
        if L > k:
            self.keys(counts, ['anchors', 'paths', 'labels', 'occurrences'] + selection_counts
                      + (['supported_paths'] if supported else []), path + '.counts')
            anchors = counts['anchors']
            self.count(anchors, 'anchors', path + '.counts.anchors',
                       extra=['by_strand' if strand_stated else 'by_orientation'])
            self.by_orientation(anchors, e, strand_stated, 'anchors', path + '.counts.anchors')
            if paths:
                self.paths_count(counts['paths'], anchors, e, strand_stated, a['limits'],
                                 path + '.counts.paths')
            elif supported:
                self.supported_paths_count(counts, anchors, e, strand_stated, request, a,
                                           predicate, path + '.counts')
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
            self.keys(counts, ['contexts', 'labels', 'occurrences'] + selection_counts,
                      path + '.counts')
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
            if labelled and f == 'labels' and (paths or supported):
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
                  + (['annotation_rows', 'annotation_units', 'memory_bytes']
                     if labelled or supported else [])
                  + (['annotation_rows_distinct'] if labelled and not supported else [])
                  # the supported-path search (SPEC §20.6): its reads in every mode
                  + (['anchor_rows', 'row_cache_hits', 'row_cache_evictions'] if supported
                     else [])
                  + (['mirror_rows'] if supported and predicate else [])
                  + (['extension_edges', 'extension_anchors', 'extension_branches']
                     if paths or supported else [])
                  + (['verification_steps'] if labelled and paths else [])
                  # a predicate: the selection's work, and the account (every mode)
                  + (['predicate_rows', 'predicate_units', 'predicate_lookups']
                     + ([] if labelled or supported else ['memory_bytes']) if predicate else [])
                  + (['motif_units'] if motif else []),
                  path + '.work')
        self.ok(all(is_int(v) and v >= 0 for v in w.values()), path + '.work')
        self.ok(w['steps'] >= w['ranges_visited'] + w.get('extension_edges', 0), path + '.work',
                'steps >= ranges_visited (+ extension_edges, one step each)')
        # SPEC §8.7: the counters beside the steps
        if paths or supported:
            ext = counts['paths']['extension'] if paths else counts['supported_paths']['search']
            examined = counts['paths' if paths else 'supported_paths']['candidates_examined']
            if ext in ('no_anchors', 'not_started', 'not_admitted'):
                self.ok((w['extension_edges'], w['extension_anchors'], w['extension_branches'])
                        == (0, 0, 0), path + '.work', f'an extension {ext}: it did not run')
            if ext == 'completed':
                self.ok(w['extension_anchors'] == anchors['value'], path + '.work',
                        'a completed extension began at every anchor')
                # each branching enters two candidates or more, its own
                self.ok(2 * w['extension_branches'] <= examined,
                        path + '.work', 'two candidates per branching')
            if anchors['relation'] == 'exact':
                self.ok(w['extension_anchors'] <= anchors['value'], path + '.work',
                        'at most the anchors listed')
        if supported:
            # SPEC §20.6: the anchors' whole rows and the mirror walks' are rows read; the
            # mirrors are read for "either" on a BASIC graph with one strand searched (§20.9)
            self.ok(w['anchor_rows'] <= w['annotation_rows'], path + '.work.anchor_rows',
                    'the anchors\' whole rows are rows read')
            if predicate:
                self.ok(w['mirror_rows'] <= w['annotation_rows'], path + '.work.mirror_rows',
                        'the mirror walks\' rows are rows read')
                self.ok(w['mirror_rows'] == 0
                        or (a['index']['graph_mode'] == 'basic'
                            and a['predicate']['strands'] == 'either'
                            and len(e['strands']) == 1 and not e['palindromic']),
                        path + '.work.mirror_rows', 'mirror walks read: "either", one strand')
        if labelled and not supported:
            # (with "predicate_only" the rows of the selected contexts were read by the
            # selection, work.predicate_rows: the projection's own reads are the placement's)
            self.ok(w['annotation_rows_distinct']
                    <= w['annotation_rows'] + w.get('predicate_rows', 0), path + '.work',
                    'each distinct row read is a row read')
            if paths and e['placement'] in ('none', 'none_canonical', 'not_requested'):
                self.ok(w['verification_steps'] == 0, path + '.work',
                        'no coordinates read: nothing verified')
        stop = e['stop']
        if stop is not None:
            self.keys(stop, SCHEMA['stop'], path + '.stop')
            self.one_of(stop['phase'], STOP_PHASES, path + '.stop.phase')
            self.one_of(stop['reason'], STOP_REASONS, path + '.stop.reason')
            self.ok(total['relation'] != 'exact'
                    or stop['phase'] in ('extraction', 'extension', 'selection')
                    + ANNOTATION_PHASES,
                    path + '.stop', 'a stop in discovery leaves no exact count')
            self.ok(stop['phase'] not in ANNOTATION_PHASES or labelled
                    or ((predicate or supported) and stop['phase'] == 'output'),
                    path + '.stop.phase', 'an annotation phase without labels "all"')
            # SPEC §19.8: the selection's phase and reasons with a predicate only
            self.ok(stop['phase'] != 'selection' or predicate, path + '.stop.phase',
                    'a selection without a predicate')
            self.ok(stop['reason'] not in ('max_predicate_work', 'max_predicate_contexts')
                    or predicate, path + '.stop.reason', 'a predicate\'s reason without one')
            self.ok(stop['reason'] != 'max_predicate_work' or stop['phase'] == 'selection',
                    path + '.stop', 'max_predicate_work stops the selection')
            self.ok(stop['reason'] != 'max_predicate_contexts'
                    or stop['phase'] in ('discovery', 'extension'), path + '.stop',
                    'max_predicate_contexts stops the raw discovery')
            self.ok(stop['phase'] != 'extension' or paths or supported, path + '.stop.phase',
                    'an extension without long_search "paths" or "supported_paths"')
            self.ok((stop['reason'] == 'max_paths') <= (stop['phase'] == 'extension'),
                    path + '.stop.reason', 'max_paths stops the extension only')
            # SPEC §20.5: the supported-path search's reads stop the extension (and a
            # predicate's mirror reads its selection); the annotation budgets elsewhere are the
            # labels' phases
            self.ok(stop['reason'] not in ('max_annotation_work', 'max_memory')
                    or stop['phase'] not in ('extension', 'selection') or supported
                    or (predicate and stop['phase'] == 'selection'
                        and stop['reason'] == 'max_memory'),
                    path + '.stop', 'an annotation budget stopping a search that reads none')
            if stop['phase'] == 'extension':
                self.ok(total['relation'] == 'exact'
                        and counts['paths']['relation'] == 'at_least', path + '.stop',
                        'a stop in the extension: the anchors exact, the paths at_least')
                if supported:
                    self.ok(counts['supported_paths']['relation'] == 'at_least'
                            and counts['supported_paths']['search'] == 'stopped',
                            path + '.stop', 'a stop in the search: the supported paths at_least')
            if stop['reason'] == 'max_steps':
                self.ok(w['steps'] <= a['limits']['max_steps'], path + '.work.steps')
            if stop['reason'] in ('max_contexts', 'max_anchors', 'max_predicate_contexts') \
                    and self.fraction is not None and stop['phase'] != 'selection':
                # stop_at_threshold without the mask compares the running upper bound: the
                # count left is above the threshold, or the note says that it may not be (the
                # selection's max_contexts compares the selected count)
                self.ok(total['value'] > a['limits'][stop['reason']]
                        or 'threshold_upper_bound' in e['notes'], path + '.notes',
                        'a threshold stop below the threshold without threshold_upper_bound')
        else:
            # without the mask a completed discovery leaves bounds where source dummies may be
            # among the candidates (SPEC §7.4)
            self.ok(total['relation'] == 'exact'
                    or (self.fraction is not None and total['relation'] == 'bounds'),
                    path + '.counts', 'a count without a stop is exact (bounds only on a graph '
                    'without its mask)')
            if paths:
                self.ok(counts['paths']['relation'] == 'exact'
                        or counts['paths']['extension'] == 'not_admitted', path + '.counts.paths',
                        'paths without a stop: exact, or not admitted')
            if supported:
                sp = counts['supported_paths']
                self.ok(sp['search'] in ('completed', 'no_anchors', 'not_admitted'),
                        path + '.counts.supported_paths.search', 'no stop: the search ended')
                self.ok(sp['relation'] == 'exact' or sp['search'] == 'not_admitted'
                        or sp.get('branches_pruned_by_predicate', 0) > 0,
                        path + '.counts.supported_paths', 'supported paths without a stop: '
                        'exact, but where a predicate pruned, or not admitted')
        self.one_of(e['determinism'], ('full', 'time_limited'), path + '.determinism')
        # SPEC §7.6/§7.9: time_limited iff the clock touched the entry, in its stop or, in
        # partial, only in its cut (stop keeps the first stop: a pattern stopped by max_steps or
        # its threshold whose release met the work time, or a later pattern after a sticky
        # max_steps stop whose empty release met it). With labels "all", a time stop of the
        # reads after a step stop of the engine states the engine's stop, so only that
        # direction is checked there.
        clocked = (stop is not None and stop['reason'] == 'time') \
            or (e.get('cut') or {}).get('reason') == 'time'
        if clocked:
            self.ok(e['determinism'] == 'time_limited', path + '.determinism',
                    'a time stop or a time cut is not deterministic')
        # SPEC §7.8: the engine's diagnostic runs on a search that completed or stopped at a
        # threshold (the engine's, or on supported paths a predicate's raised through its
        # selection); a stop of a later phase comes after it and leaves the note
        threshold_stop = stop is not None and (
            (stop['phase'] in ('discovery', 'extension')
             and stop['reason'] in ('max_contexts', 'max_anchors', 'max_paths',
                                    'max_predicate_contexts'))
            or stop['phase'] in ('selection', 'output') + ANNOTATION_PHASES)
        if e['determinism'] == 'time_limited':
            # (SPEC §7.8: the low-complexity diagnostic of a pattern of more than 191 bases, cut
            # by the work time, leaves its note out and states time_limited, its stop unchanged)
            diagnosis_cut = (stop is None or threshold_stop) and L > LOW_COMPLEXITY_UNCUT \
                and 'low_complexity_pattern' not in e['notes']
            self.ok(clocked or ((labelled or predicate or supported) and stop is not None)
                    or diagnosis_cut, path + '.stop',
                    'time_limited without a time stop or a time cut')
        self.ok(isinstance(e['notes'], list) and all(n in NOTES for n in e['notes']),
                path + '.notes', repr(e['notes']))
        self.ok(e['notes'] == [n for n in NOTES if n in e['notes']], path + '.notes', 'order')
        # SPEC §7.8: never beside a budget stop of the search
        self.ok('low_complexity_pattern' not in e['notes'] or stop is None or threshold_stop,
                path + '.notes', 'low_complexity_pattern beside a budget stop')
        self.ok(('strand_unknown_canonical' in e['notes']) is (not strand_stated),
                path + '.notes', 'strand_unknown_canonical exactly where no strand is known')
        self.ok(('paths_later_increment' in e['notes']) is (L > k and not paths
                                                            and not supported),
                path + '.notes', 'paths_later_increment exactly for L > k without paths')
        # (a predicate answer never says annotation_not_read: its selection read rows; it says
        # projection_not_read for a projection named and not built, SPEC §19.10; so does a
        # pattern searched for its supported paths, whose search read rows, SPEC §20.6)
        self.ok(('annotation_not_read' in e['notes']) is (named and not labelled
                                                          and not predicate and not supported),
                path + '.notes', 'annotation_not_read exactly where labels were named, not read')
        if supported:
            self.ok(('projection_not_read' in e['notes'])
                    is (supported_projection_named and not labelled), path + '.notes',
                    'projection_not_read exactly where a projection of the supported paths was '
                    'named and not built')
        else:
            self.ok(('projection_not_read' in e['notes'])
                    is (predicate and projection_named and not labelled), path + '.notes',
                    'projection_not_read exactly where a predicate request named a projection '
                    'it did not build')
        if predicate:
            self.ok(('predicate_constant' in e['notes'])
                    is (e['selection']['pass'] == 'constant'), path + '.notes',
                    'predicate_constant exactly for a constant normal form')
        else:
            self.ok('predicate_constant' not in e['notes'], path + '.notes')
        # graphs without their mask, and the stop '*'
        self.ok(('estimate_sampled_dummy_fraction' in e['notes']) is self.estimated,
                path + '.notes', 'estimate_sampled_dummy_fraction exactly where a count states an '
                'estimate')
        self.ok(('no_stop_codon' in e['notes']) is (not pk.has_instances()), path + '.notes',
                'no_stop_codon exactly for a peptide whose * has no codon in its table')
        # SPEC §7.4: a pattern whose discovery completed with at most
        # caps.max_checked_entries unchecked candidates had each of them tested and is exact,
        # so a bounds total stated without a stop has more of them: on a BASIC or CANONICAL
        # graph they number its upper - lower (a wrapped PRIMARY graph counts an entry in both
        # orientations; partial's release raises the lower bounds to what it listed)
        if self.checked_limit is not None and e['stop'] is None \
                and self.graph_mode != 'primary' and mode != 'partial' \
                and total['relation'] == 'bounds':
            self.ok(total['upper'] - total['lower'] > self.checked_limit, path + '.counts',
                    f'bounds with {total["upper"] - total["lower"]} unchecked candidates, at most '
                    f'max_checked_entries ({self.checked_limit}): they are checked, exact')
        if not pk.has_instances():
            if L <= k:
                self.ok(total['relation'] == 'exact' and total['value'] == 0, path + '.counts',
                        'no instance: exact 0')
            elif paths and counts['paths']['relation'] in ('exact', 'at_least'):
                self.ok(counts['paths']['value'] == 0, path + '.counts.paths',
                        'no instance: no path')
        if 'threshold_upper_bound' in e['notes']:
            # a threshold decided on the upper bound (no mask) against the request while the
            # lower bound did not cross it: the decision is stated (SPEC §7.5)
            ext = counts.get('paths', {}).get('extension')
            self.ok(self.fraction is not None, path + '.notes',
                    'threshold_upper_bound only on a graph without its mask')
            self.ok((e.get('withheld') or {}).get('reason') in ('count_above_threshold',
                                                                 'threshold_crossed',
                                                                 'anchors_above_threshold',
                                                                 'predicate_above_threshold')
                    or (stop or {}).get('reason') in ('max_contexts', 'max_anchors',
                                                      'max_predicate_contexts')
                    or ext == 'not_admitted'
                    or (e.get('selection') or {}).get('pass') == 'not_admitted',
                    path + '.notes', 'threshold_upper_bound without a threshold decision')
        if labelled:
            self.ok(('annotation_unbudgeted' in e['notes'])
                    is (e['annotation'] == 'unbudgeted'), path + '.notes')
            self.ok(('record_bounds_unknown' in e['notes']) is (e['placement'] == 'global'),
                    path + '.notes')
            # the labels of paths where no coordinate is read; of supported paths, searched at
            # the label level where no occurrence is listed (SPEC §20.6)
            if supported:
                expected = counts['supported_paths']['level'] == 'label_intersection' \
                    and e['placement'] != 'global'
            else:
                expected = paths and e['placement'] in ('none', 'none_canonical',
                                                        'not_requested')
            self.ok(('label_intersection_only' in e['notes']) is expected, path + '.notes',
                    'label_intersection_only exactly for paths placed nowhere')
        else:
            self.ok(not {'record_bounds_unknown', 'label_intersection_only'} & set(e['notes']),
                    path + '.notes')
            # a selection (or a supported-path search) that read rows of an unbudgeted
            # annotation says so
            unbudgeted = 'annotation_unbudgeted' in e['notes']
            read = (predicate and e['work']['predicate_rows'] > 0) \
                or (supported and e['work']['annotation_rows'] > 0)
            self.ok(not unbudgeted or read, path + '.notes',
                    'annotation_unbudgeted where nothing was read')
            if read and self.capabilities is not None:
                self.ok(unbudgeted is (self.capabilities['annotation'] == 'unbudgeted'),
                        path + '.notes', 'annotation_unbudgeted exactly on an unbudgeted index')
        self.timing(e['timing'], path + '.timing', labelled, paths, predicate,
                    supported=supported)
        if predicate:
            self.selection(e, total, L, k, request, a, path)
            self.rows_refused(e, path)
        elif supported:
            self.rows_refused(e, path)
        if motif:
            self.motif(e, total, L, k, a, path)

        if mode == 'count':
            self.ok(e['retrieval_complete'] is False, path + '.retrieval_complete',
                    'a count returns no context')
            return
        self.retrieval(e, pk, total, request, a, path, paths, supported)

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
            # P[L - k, L)'s; a peptide's windows may cut a codon (exact bits, SPEC §12.2)
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

    def retrieval(self, e, pk, total, request, a, path, paths=False, supported=False):
        k = a['index']['k']
        strand_stated = a['index']['strand_stated']
        mode = a['mode']
        labelled = isinstance(a['output'], dict) and a['output']['labels'] in ('all',
                                                                               'predicate_only')
        predicate = 'predicate' in request
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
        if not predicate:
            # SPEC §19.8, §20.9: the selection's reasons come with a predicate only
            self.ok((e['withheld'] or {}).get('reason') not in (
                        'selected_above_threshold', 'predicate_above_threshold',
                        'predicate_budget')
                    and (e['cut'] or {}).get('reason') not in ('max_predicate_contexts',
                                                               'max_predicate_work'),
                    path, 'a predicate\'s withheld or cut without a predicate')
        if e['cut'] is not None:
            self.one_of(e['cut']['reason'], CUT, path + '.cut.reason')
            self.ok(mode == 'partial' and not e['retrieval_complete'], path + '.cut',
                    'only partial cuts, and a cut list is not complete')
            self.ok(e['cut']['reason'] != 'max_contexts' or L <= k, path + '.cut',
                    'max_contexts cuts contexts')
            self.ok(e['cut']['reason'] != 'max_paths'
                    or ((paths or supported) and e['returned'] == limits['max_paths']),
                    path + '.cut', 'max_paths: the first max_paths paths')
            # SPEC §20.5: the supported-path search's own work budget
            self.ok(e['cut']['reason'] != 'max_annotation_work' or supported, path + '.cut',
                    'max_annotation_work cuts the supported paths')
        if e['retrieval_complete'] and predicate:
            # SPEC §19.11: every selected context returned (a constant false selects none,
            # whatever stopped the raw search)
            sel = e['counts']['selected']
            self.ok(sel['relation'] == 'exact' and e['returned'] == sel['value']
                    and e['withheld'] is None and e['cut'] is None, path,
                    'complete with a predicate: selected exact, every selected context returned')
            self.ok(e['selection']['pass'] in ('completed', 'constant', 'not_started'), path,
                    'complete: the selection completed (or needed no pass)')
            if e['selection']['pass'] != 'constant' or self.normal_form is True:
                self.ok(total['relation'] == 'exact' and e['stop'] is None, path,
                        'complete: exact, no stop')
        elif e['retrieval_complete']:
            self.ok(total['relation'] == 'exact' and e['stop'] is None
                    and e['withheld'] is None and e['cut'] is None, path,
                    'complete: exact, no stop, nothing withheld or cut')
            if paths:
                c = e['counts']['paths']
                self.ok(c['relation'] == 'exact' and e['returned'] == c['value'], path,
                        'complete: every path returned')
            elif supported:
                c = e['counts']['supported_paths']
                self.ok(c['relation'] == 'exact' and e['returned'] == c['value'], path,
                        'complete: every supported path returned')
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
        # the admissions (SPEC §7.5; without the mask they compare the upper bound, and
        # threshold_upper_bound states a decision its lower bound did not cross)
        reason = (e['withheld'] or {}).get('reason')
        noted = 'threshold_upper_bound' in e['notes']
        admission = {'count_above_threshold': 'max_contexts',
                     'anchors_above_threshold': 'max_anchors'}.get(reason)
        if admission and (L <= k or admission == 'max_anchors'):
            self.ok(self.above(total, limits[admission]), path + '.withheld',
                    f'the exact count (or its upper bound) above {admission}')
            if total['relation'] == 'bounds':
                self.ok(noted is (total['lower'] <= limits[admission]), path + '.notes',
                        f'threshold_upper_bound exactly where the lower bound is within '
                        f'{admission}')
        if L > k and not paths and not supported:
            self.ok(results == [], path + '.results', 'no result without long_search "paths"')
            if not e['retrieval_complete']:
                self.ok(e['withheld'] == {'reason': 'paths_later_increment'}, path + '.withheld')
        if paths:
            self.ok(e['withheld'] != {'reason': 'paths_later_increment'}, path + '.withheld',
                    'paths were asked for')
            c = e['counts']['paths']
            reason = (e['withheld'] or {}).get('reason')
            if reason == 'anchors_above_threshold':
                self.ok(self.above(total, limits['max_anchors'])
                        and c['extension'] == 'not_admitted', path + '.withheld',
                        'the exact anchors (or their upper bound) above max_anchors, the '
                        'extension not admitted')
            elif reason == 'count_above_threshold':
                self.ok(c['relation'] == 'exact' and c['value'] > limits['max_paths'],
                        path + '.withheld', 'the exact paths above max_paths')
            elif reason == 'threshold_crossed':
                self.ok(e['stop'] is not None
                        and e['stop']['reason'] in ('max_anchors', 'max_paths'),
                        path + '.withheld', 'a stop_at_threshold stop')
        if supported:
            # SPEC §20.5: the release of the supported paths and its reasons
            self.ok(e['withheld'] != {'reason': 'paths_later_increment'}, path + '.withheld',
                    'supported paths were asked for')
            c = e['counts']['supported_paths']
            stop = e['stop'] or {}
            reason = (e['withheld'] or {}).get('reason')
            if reason == 'anchors_above_threshold':
                self.ok(self.above(total, limits['max_anchors'])
                        and c['search'] == 'not_admitted', path + '.withheld',
                        'the exact anchors (or their upper bound) above max_anchors, the '
                        'search not admitted')
            elif reason == 'count_above_threshold':
                self.ok(not predicate and c['relation'] == 'exact'
                        and c['value'] > limits['max_paths'], path + '.withheld',
                        'the exact supported paths above max_paths')
            elif reason == 'threshold_crossed':
                self.ok(stop.get('reason') in ('max_anchors', 'max_paths',
                                               'max_predicate_contexts'),
                        path + '.withheld', 'a stop_at_threshold stop')
            elif reason == 'annotation_budget':
                self.ok(stop.get('reason') in ('max_annotation_work', 'max_memory')
                        or e['rows_refused'], path + '.withheld',
                        'the search\'s (or the mirrors\') reads stopped')
            elif reason == 'output_budget':
                self.ok(stop.get('reason') == 'max_memory', path + '.withheld',
                        'the list did not fit the account')
            if (e['cut'] or {}).get('reason') == 'max_annotation_work':
                self.ok(stop.get('reason') == 'max_annotation_work', path + '.cut',
                        'a stop at max_annotation_work')
            if (e['cut'] or {}).get('reason') == 'max_memory':
                # the account ended the list (stop {output, max_memory}) or the search's or the
                # selection's reads; the first stop is the one stated
                self.ok(e['stop'] is not None, path + '.cut',
                        'a max_memory cut states its stop')
        if mode == 'partial':
            self.ok(e['returned'] <= (limits['max_paths'] if paths or supported
                                      else limits['max_contexts']),
                    path + '.returned', 'the cap')
        if predicate:
            self.predicate_retrieval(e, total, L, k, a, path, supported)

        if paths or supported:
            self.path_results(e, pk, request, a, path, supported)
        else:
            self.context_results(e, pk, a, path)
        if labelled:
            self.labels(e, total, L, a, path, paths or supported, request)

    def predicate_retrieval(self, e, total, L, k, a, path, supported=False):
        """What a predicate's withheld, cut and results mean (SPEC §19.8, §19.10; on supported
        paths §20.9: the walks held for "either" are the raw ones, max_paths the selected
        list's threshold)."""
        limits = a['limits']
        s = e['selection']
        sel = e['counts']['selected']
        stop = e['stop'] or {}
        reason = (e['withheld'] or {}).get('reason')
        if reason == 'predicate_above_threshold':
            if supported:
                held = e['counts']['supported_paths']
                self.ok(s['pass'] in ('not_admitted', 'not_started')
                        and held['relation'] != 'unknown'
                        and held['value'] > limits['max_predicate_contexts'], path + '.withheld',
                        'more supported walks to hold than max_predicate_contexts')
            else:
                self.ok(s['pass'] == 'not_admitted'
                        and self.above(total, limits['max_predicate_contexts']),
                        path + '.withheld', 'the raw count (or its upper bound) above '
                        'max_predicate_contexts, not admitted')
        if reason == 'selected_above_threshold':
            cap = limits['max_paths'] if supported else limits['max_contexts']
            self.ok(sel['relation'] == 'exact' and sel['value'] > cap,
                    path + '.withheld', 'the exact selected count above max_contexts '
                    '(max_paths)')
        if reason == 'predicate_budget':
            self.ok(s['pass'] in ('stopped', 'not_started')
                    and (stop.get('reason') in ('max_predicate_work', 'max_memory')
                         or e['rows_refused']), path + '.withheld',
                    'the selection stopped by its work or the account')
        if reason == 'threshold_crossed':
            self.ok(stop.get('reason') in ('max_contexts', 'max_predicate_contexts')
                    + (('max_paths', 'max_anchors') if supported else ()),
                    path + '.withheld', 'a stop_at_threshold stop of the pass or the raw search')
        cut = (e['cut'] or {}).get('reason')
        if cut == 'max_predicate_work':
            self.ok(stop.get('reason') == 'max_predicate_work', path + '.cut')
        if cut == 'max_predicate_contexts':
            # the raw release's cut applies only when neither the engine nor the pass stopped
            # (their cuts come first): every admitted raw context was decided
            self.ok(stop.get('reason') in (None, 'max_predicate_contexts')
                    and e['counts']['tested']['relation'] == 'exact'
                    and e['counts']['tested']['value'] <= limits['max_predicate_contexts'],
                    path + '.cut', 'the raw release held the first max_predicate_contexts')
            # (supported paths: a walk held whose mirror is unknown may stay undecided)
            if not stop and not supported:
                self.ok(e['counts']['tested']['value'] == limits['max_predicate_contexts'],
                        path + '.cut', 'the raw release held the first max_predicate_contexts')
        if s['pass'] == 'not_admitted' and a['mode'] != 'count':
            self.ok(reason == 'predicate_above_threshold' or e['returned'] == 0, path,
                    'not admitted: nothing listed')
        if L <= k or supported:
            # a supported path's own support is its labels at the level searched (SPEC §20.9)
            level = e['counts']['supported_paths']['level'] if supported else None
            labels = isinstance(a['output'], dict) and a['output']['labels'] in (
                'all', 'predicate_only')
            if labels and e['results'] and self.normal_form is not None:
                self.selection_labels(e, path, a, level)
            if a['output'] == {'labels': 'predicate_only',
                               'occurrences': a['output'].get('occurrences')} \
                    and isinstance(self.normal_form, dict):
                names = set(predicate_names(self.normal_form))
                for i, r in enumerate(e['results']):
                    for label in r['labels'] or []:
                        self.ok(label['column'] in names, f'{path}.results[{i}].labels',
                                'predicate_only: the predicate\'s labels only')
                        if level == 'record_verified' and label['support'] != 'record_verified':
                            continue
                        self.ok(label['column'] in r['selection_labels'],
                                f'{path}.results[{i}].labels',
                                'a label of its own row is in the set it was evaluated on')

    def path_results(self, e, pk, request, a, path, supported=False):
        """The paths of long_search "paths" (SPEC §12.1): each with its L bases (sequence =
        instance, an instance of the oriented pattern), its anchor's k bases, offset 0, its node
        path (n = L - k + 1 node ids, the row of each k-mer); never kmer; in the answer order
        (anchor node, orientation, sequence), each once; per strand within the counts."""
        k = a['index']['k']
        strand_stated = a['index']['strand_stated']
        labelled = isinstance(a['output'], dict) and a['output']['labels'] in ('all',
                                                                               'predicate_only')
        # a selected supported path (SPEC §20.9): why it was selected, under a projection
        predicate = 'predicate' in request
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
                      + (['labels_excluded_unverified'] if labelled and verified_only else [])
                      + (['selection_labels', 'selection_strands'] if labelled and predicate
                         else []), rp)
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
        c = e['counts']['supported_paths' if supported else 'paths']
        by = c['by_strand' if strand_stated else 'by_orientation']
        got = {}
        for r in e['results']:
            part = STRAND_KEYS[r[key]] if strand_stated else r[key]
            got[part] = got.get(part, 0) + 1
        for part, count in by.items():
            m = got.get(part, 0)
            if e['retrieval_complete'] and not predicate:
                self.ok(count['relation'] == 'exact' and m == count['value'],
                        f'{path}.counts.paths.{part}', f'{m} results, the count {count["value"]}')
            elif count['relation'] == 'exact':
                self.ok(m <= count['value'], f'{path}.counts.paths.{part}',
                        f'{m} results above the exact count {count["value"]}')
        self.ok(set(got) <= set(by), path + '.counts.paths', 'a result outside the parts counted')

    def context_results(self, e, pk, a, path):
        k = a['index']['k']
        strand_stated = a['index']['strand_stated']
        labelled = isinstance(a['output'], dict) and a['output']['labels'] in ('all',
                                                                               'predicate_only')
        # the results of a predicate request are its selected contexts, with their
        # selection_labels when the projection reads labels
        predicate = 'predicate' in a
        L = pk.length()
        results = e['results']
        key = 'strand' if strand_stated else 'orientation'
        allowed = set(e['strands'])
        order = []
        for i, r in enumerate(results):
            rp = f'{path}.results[{i}]'
            self.keys(r, ['kmer', 'instance', 'offset', key, 'node', 'row']
                      + (['support', 'labels_status', 'labels_total', 'labels'] if labelled
                         else [])
                      + (['selection_labels', 'selection_strands'] if labelled and predicate
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
        # the results and the counts beside them
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
                    if e['retrieval_complete'] and not predicate:
                        # every context returned: the results per part are the part's count
                        self.ok(count['relation'] == 'exact' and n == count['value'],
                                f'{path}.counts.contexts.{name}.{part}',
                                f'{n} results, the count {count["value"]}')
                    elif count['relation'] == 'exact':
                        self.ok(n <= count['value'], f'{path}.counts.contexts.{name}.{part}',
                                f'{n} results above the exact count {count["value"]}')
                    elif count['relation'] == 'bounds' and self.fraction is not None:
                        # without the mask a released context is a k-mer among the upper
                        # bound's candidates, and the release raised the lower bound to it
                        self.ok(n <= count['lower'], f'{path}.counts.contexts.{name}.{part}',
                                f'{n} results above the lower bound {count["lower"]}')
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
        self.rows_refused(e, path)
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
            # partial whose memory account could not hold by_label (SPEC §14.4): no label of
            # the pattern is built, every context read answers output_budget, the answer is
            # incomplete and states a stop (the output's, or an earlier one: the first stop
            # wins); the counts keep their own relations
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
        """The labels of paths (SPEC §12.1): by_label per label with its paths, the paths one
        record verifies and its occurrences; each path's labels, each with its support
        (record_verified: an occurrence of the whole path in one record, placed;
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

    def block(self, b, path, server_wide=False):
        """The full block (SPEC §10.2). |server_wide|: GET /capabilities of a multi-graph
        server (SPEC §24), the contract and the caps without a graph (its graph fields null,
        as while a single index loads), available whether a pair is served"""
        self.keys(b, SCHEMA['capabilities'], path)
        if server_wide:
            self.ok(b['available'] in (True, False), path + '.available', 'known')
            for f in GRAPH_FIELDS:
                self.ok(b[f] is None, f'{path}.{f}', 'per pair on a multi-graph server')
            self.ok(b['predicate']['access'] is None, path + '.predicate.access')
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
        # the paths are opt-in: an omitted long_search is "anchors", and long_patterns keeps
        # describing that answer
        self.ok(b['default_long_search'] == 'anchors' and 'anchors' in b['long_search'],
                path + '.default_long_search', 'version 1: an omitted long_search is "anchors"')
        self.ok(b['long_patterns'] == 'anchors_counted', path + '.long_patterns')
        # peptides: the residues, the genetic codes (NCBI ids) and the default
        self.ok(isinstance(b['protein_residues'], list)
                and all(isinstance(r, str) and len(r) == 1 and r in PROTEIN_RESIDUES
                        for r in b['protein_residues']), path + '.protein_residues')
        self.ok(isinstance(b['genetic_codes'], list) and b['genetic_codes']
                and all(c in GENETIC_CODES for c in b['genetic_codes']), path + '.genetic_codes',
                'NCBI translation table ids')
        self.ok(b['default_genetic_code'] in b['genetic_codes'], path + '.default_genetic_code')
        # SPEC §10.2: the two prose fields are references to the SPEC, in printable ASCII (the
        # server escapes any other byte as \uXXXX), and the delivery rates of the time kept
        # back for the answer (§7.6) are numbers
        for f in ('protein_rule', 'caps_rule'):
            self.ok(isinstance(b[f], str) and 'SPEC-pattern-search.md section' in b[f]
                    and all(' ' <= c <= '~' for c in b[f]), f'{path}.{f}',
                    'a reference to the SPEC, printable ASCII')
        self.keys(b['delivery_mbps'], SCHEMA['delivery_mbps'], path + '.delivery_mbps')
        self.ok(all(is_num(v) and v > 0 for v in b['delivery_mbps'].values()),
                path + '.delivery_mbps', 'positive rates in MB/s')
        self.ok(('protein' in b['kinds']) is bool(b['protein_residues']), path + '.kinds',
                'protein served exactly with its residues')
        self.ok(b['modes'] == list(MODES), path + '.modes')
        self.ok('none' in b['projections'], path + '.projections', '"none" is served everywhere')
        self.ok(b['default_mode'] in b['modes'], path + '.default_mode')
        self.ok(b['default_projection'] in b['projections'], path + '.default_projection')
        # contract version 1: an omitted output.labels means "none" on every server stating
        # it (a "may become all" would change the meaning of an unchanged request under the
        # same version)
        self.ok(b['default_projection'] == 'none', path + '.default_projection',
                'version 1: an omitted output.labels is "none"')
        # the caps are classified in SPEC §4.5 (test_caps_rule_section_classifies_every_cap)
        self.ok(b['caps_rule'] == CAPS_RULE, path + '.caps_rule', f'a reference to {CAPS_RULE}')
        self.ok(not set(b['projections']) & set(b['projections_later_increment']),
                path + '.projections_later_increment')
        self.ok(not set(b['kinds']) & set(b['kinds_later_increment']),
                path + '.kinds_later_increment')
        self.ok(b['default_scope'] in ('suffix', 'any_offset'), path + '.default_scope')
        self.ok(b['default_strands'] in b['strands'], path + '.default_strands')
        self.keys(b['scopes_by_graph_mode'], GRAPH_MODES, path + '.scopes_by_graph_mode')
        self.keys(b['caps'], CAPS, path + '.caps')
        self.ok(all(is_num(v) for v in b['caps'].values()), path + '.caps')
        # SPEC §19.12: "predicate_only" served with the predicate object (the operators, the
        # strands, the access: null while the annotation is not described)
        self.keys(b['predicate'], SCHEMA['capabilities_predicate'], path + '.predicate')
        self.ok(b['predicate']['operators'] == list(PREDICATE_OPERATORS),
                path + '.predicate.operators')
        self.ok(b['predicate']['strands'] == list(PREDICATE_STRANDS), path + '.predicate.strands')
        # SPEC §25.7: predicate_scope's values
        self.ok(b['predicate']['scopes'] == list(PREDICATE_SCOPE_VALUES),
                path + '.predicate.scopes')
        self.ok(b['predicate']['access'] is None or b['predicate']['access'] in ACCESS,
                path + '.predicate.access')
        self.ok(('predicate_only' in b['projections'])
                is ('predicate_only' not in b['projections_later_increment']),
                path + '.projections', 'predicate_only served or announced')
        self.ok(1 <= b['caps']['max_predicate_work'] and b['caps']['max_predicate_labels']
                <= 1_000_000, path + '.caps', 'the predicate caps in their ranges')
        self.ok(is_num(b['default_time_budget_ms'])
                and b['finalize_reserve_ms'] < b['default_time_budget_ms']
                <= b['caps']['time_budget_ms'], path + '.default_time_budget_ms')
        # the defaults of max_anchors and max_labels (SPEC §4.5): integers at most their caps
        # (the server refuses a larger default at start-up), not among the caps
        for cap in ('max_anchors', 'max_labels'):
            self.ok(is_int(b['default_' + cap]) and 0 <= b['default_' + cap] <= b['caps'][cap]
                    and 'default_' + cap not in b['caps'], f'{path}.default_{cap}',
                    'at most caps.' + cap)
        self.ok(isinstance(b['caps_rule'], str), path + '.caps_rule')
        # in_ram is accepted, as /search accepts it (SPEC §24): the index is loaded for the
        # request where the server runs on mmap
        self.ok(b['resident_only'] is False and b['in_ram'] == 'accepted',
                path + '.resident_only', 'in_ram accepted')
        if server_wide:
            # the graph fields are the pairs' (graph_summary, the pairs' own blocks)
            return
        # SPEC §10.2: how an available graph counts, and the dummy fraction its estimates rest
        # on; null when the graph is not served (or loading)
        if b['available'] is True:
            self.one_of(b['counting'], COUNTINGS, path + '.counting')
            self.ok(b['counting'] == ('upper_bound' if b['mask'] == 'absent' else 'exact'),
                    path + '.counting', 'exact with a mask, upper_bound without one')
            if b['counting'] == 'upper_bound':
                self.dummy_fraction(b['dummy_fraction'], path + '.dummy_fraction')
            else:
                self.ok(b['dummy_fraction'] is None, path + '.dummy_fraction',
                        'null with the mask')
        else:
            for f in ('counting', 'dummy_fraction'):
                self.ok(b[f] is None, f'{path}.{f}', 'null on a graph not served')
        if b['available'] is None:
            for f in ('graph_mode', 'k', 'alphabet', 'strand_stated', 'mask', 'scopes',
                      'placement', 'support', 'annotation'):
                self.ok(b[f] is None, f'{path}.{f}', 'unknown while the index loads')
            self.ok(b['predicate']['access'] is None, path + '.predicate.access')
            return
        self.ok(is_int(b['k']), path + '.k')
        if b['graph_mode'] is None:
            # SPEC §10.2: a graph the engine does not recognise: not available, and only k is
            # set
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
        # mask_required is retired: a block that states counting is never refused for a
        # missing mask
        self.ok(b['unavailable_reason'] not in RETIRED, path + '.unavailable_reason',
                'a retired reason (mask_required: a graph without its mask is served)')
        # a mask that marks a W = $ edge valid is a mask read from its file (one built at load
        # is not checked, SPEC §6)
        if b['unavailable_reason'] == 'mask_invalid':
            self.ok(b['mask'] == 'file', path + '.mask', 'mask_invalid: a mask file')
        # the alphabet before the mask (pattern.cpp route_support): on another alphabet than the
        # served one the alphabet's reason is given, mask or none -- a DNA5 graph without a mask
        # is alphabet_untested
        if b['alphabet'] == SERVED_ALPHABET:
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
        # a budget-aware annotation is read by rows; an unbudgeted one with direct access by
        # single cells (at most 16 labels)
        if b['annotation'] == 'budgeted':
            self.ok(b['predicate']['access'] == 'rows', path + '.predicate.access')
        if b['graph_mode'] != 'basic' and b['placement'] is not None:
            self.ok(b['placement'] == 'none_canonical', path + '.placement')


class Raising:
    """The Checker's test case for the mutation tests: a rule that fails raises."""

    def fail(self, message):
        raise AssertionError(message)


class TestPatternFixtures(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.index = load('index.json')
        cls.fixtures = cls.index['fixtures']
        cls.bodies = {name: (load(name, 'request.json'), load(name, 'answer.json'))
                      for name in cls.fixtures}
        # a multi-graph server's answers (SPEC §24) as stored; the checks of every answer see
        # the answer of the fixture server's one answered pair, the pair's tags and outcome
        # taken off, which is the single-graph answer (test_multi_graph_answers checks the
        # envelope, the refused entries among its rules)
        cls.multi_answers = {}
        for name, f in cls.fixtures.items():
            request, answer = cls.bodies[name]
            if f['method'] == 'POST' and f['status'] == 200 and 'answers' in answer:
                cls.multi_answers[name] = answer
                answered = [x for x in answer['answers'] if x.get('outcome') == 'answered']
                assert len(answered) == 1, name
                pair = {k: v for k, v in answered[0].items() if k not in SCHEMA['answer_pair']}
                cls.bodies[name] = (request, pair)

    def full_block(self, name):
        """The full pattern block of GET fixture |name| (SPEC §23): the document of
        /pattern/capabilities, the `pattern` member of /capabilities, that of
        /traverse/capabilities without `details`; None for a refusal."""
        f = self.fixtures[name]
        doc = self.bodies[name][1]
        if f['status'] != 200:
            return None
        if f['path'].split('?')[0] == '/pattern/capabilities':
            # a multi-graph server's names its pair (SPEC §24)
            return {k: v for k, v in doc.items() if k not in SCHEMA['capabilities_pair']}
        return {k: v for k, v in doc['pattern'].items() if k != 'details'}

    def capabilities_of(self, server):
        """The full pattern block of the graph a fixture's answers were computed on: the first
        GET fixture of its server that describes a graph (a multi-graph server's GET
        /capabilities describes none: its pair's block is the probe's)."""
        for name, f in self.fixtures.items():
            if f['server'] == server and f['method'] == 'GET' and not f['hand_made']:
                b = self.full_block(name)
                if b is not None and b['k'] is not None:
                    return b
        return None

    def check_stored(self, name, mutate=None):
        """Checker.answer on fixture |name|'s answer, changed by |mutate| first."""
        request, answer = self.bodies[name]
        answer = copy.deepcopy(answer)
        if mutate:
            mutate(answer)
        Checker(Raising(), name).answer(answer, request,
                                        self.capabilities_of(self.fixtures[name]['server']))

    def assertRefusesMutations(self, cases):
        """Each (fixture, mutate, words): the changed answer fails a rule whose message matches
        |words| -- the rules the stored answers pass are not vacuous."""
        for name, mutate, says in cases:
            with self.subTest(fixture=name, says=says):
                with self.assertRaisesRegex(AssertionError, says):
                    self.check_stored(name, mutate)

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
        check.ok(request.get('supported_paths_level', 'best') in SUPPORTED_PATHS_LEVELS,
                 'request.supported_paths_level')
        # SPEC §20.2: the label level with record_verified under supported_paths is refused
        check.ok(not (request.get('long_search') == 'supported_paths'
                      and request.get('supported_paths_level') == 'label_intersection'
                      and request.get('require_support') == 'record_verified'),
                 'request.supported_paths_level', 'label_intersection with record_verified')
        if 'predicate' in request:
            check.predicate_form(request['predicate'], 'request.predicate')
            check.ok(request.get('long_search', 'anchors') != 'paths', 'request.long_search',
                     'a predicate with long_search "paths" is refused')
        # SPEC §25.2: predicate_scope needs a predicate
        if 'predicate_scope' in request:
            check.one_of(request['predicate_scope'], PREDICATE_SCOPE_VALUES,
                         'request.predicate_scope')
            check.ok('predicate' in request, 'request.predicate_scope', 'needs a predicate')
        check.ok(request.get('output', {}).get('labels') != 'predicate_only'
                 or 'predicate' in request, 'request.output.labels',
                 'predicate_only needs a predicate')
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
        """A GET of the three capabilities routes (SPEC §10.1, §23): the full block on
        /capabilities and /pattern/capabilities (there the document itself), the block with
        `details` on /traverse/capabilities; /pattern/capabilities with ?graph= on a
        single-graph server a 400 {error}."""
        multi = f['server'] == 'multi'
        route = f['path'].split('?')[0]
        if f['status'] != 200:
            # a single-graph server refuses ?graph= (it has one), a multi-graph server a request
            # without it (which pair?), both as /traverse/capabilities does
            selects = bool(re.search(r'[?&]graph(_path)?=', f['path']))
            check.ok(route == '/pattern/capabilities' and f['status'] == 400
                     and selects is not multi, 'status',
                     'only ?graph= on a single-graph server, or none on a multi-graph one, is '
                     'refused')
            check.keys(doc, ['error'], 'answer')
            check.ok(isinstance(doc['error'], str)
                     and ('needs ?graph=<name>' if multi else 'single graph') in doc['error'],
                     'answer.error')
            return
        if route == '/pattern/capabilities':
            if multi:
                # the pair ?graph= selected, named (SPEC §24)
                check.ok(set(SCHEMA['capabilities_pair']) <= set(doc), 'answer',
                         'the pair named')
                check.ok(re.search(r'[?&]graph=([^&]+)', f['path']).group(1) == doc['graph'],
                         'graph', 'the name ?graph= gave')
                doc = {k: v for k, v in doc.items() if k not in SCHEMA['capabilities_pair']}
            check.block(doc, 'pattern')
            return
        b = doc['pattern']
        if route == '/traverse/capabilities':
            check.ok(b.get('details') == DETAILS, 'pattern.details', DETAILS)
            check.ok(set(GATE_KEYS) <= set(b), 'pattern', 'every gate field of the block')
            b = {k: v for k, v in b.items() if k != 'details'}
        else:
            check.ok(route == '/capabilities', 'path')
            check.ok('details' not in b, 'pattern.details', 'only on /traverse/capabilities')
        check.block(b, 'pattern', server_wide=multi and route == '/capabilities')
        if route == '/capabilities':
            check.ok(doc['mode'] == ('multi' if multi else 'single'), 'mode')
            check.ok('pattern' in doc['features'], 'features')
            check.ok(doc['routes'].get('pattern') == 'POST /pattern', 'routes.pattern')
            check.ok(doc['routes'].get('pattern_capabilities')
                     == (DETAILS + '?graph={name}[&graph_path={path}]' if multi else DETAILS),
                     'routes.pattern_capabilities')
            # the largest request body the server reads (--max-request-body-mb, MiB): an
            # integer of at least 1, or null when unlimited (the fixture servers)
            check.ok(doc['max_request_body_mb'] is None
                     or (is_int(doc['max_request_body_mb']) and doc['max_request_body_mb'] >= 1),
                     'max_request_body_mb', 'an integer of at least 1, or null')
            if not multi:
                check.ok(doc['graph_summary'] is None, 'graph_summary', 'null on one graph')
                check.ok(doc['max_graphs_without_selection'] is None,
                         'max_graphs_without_selection', 'null on one graph')
                return
            # SPEC §24.1: the threshold of the selection rule (a request without `graphs`
            # queries every name of a list of at most that many), the server flag, at least 1
            check.ok(is_int(doc['max_graphs_without_selection'])
                     and doc['max_graphs_without_selection'] >= 1,
                     'max_graphs_without_selection', 'an integer of at least 1')
            # SPEC §24: every pair of the list, what the pattern search and a traversal make
            # of it, and whether their columns are disjoint
            summary = doc['graph_summary']
            check.keys(summary, SCHEMA['graph_summary'], 'graph_summary')
            check.ok(isinstance(summary['columns_disjoint'], bool)
                     and is_int(summary['shared_columns'])
                     and summary['columns_disjoint'] is (summary['shared_columns'] == 0),
                     'graph_summary', 'disjoint exactly when no column is shared')
            names = [p['graph'] for p in summary['pairs']]
            check.ok(names == sorted(names) and set(names) == set(doc['graphs']),
                     'graph_summary.pairs', 'every name of the list, in order')
            for i, p in enumerate(summary['pairs']):
                path = f'graph_summary.pairs[{i}]'
                check.keys(p, SCHEMA['graph_summary_pair'], path)
                check.keys(p['traversal'], SCHEMA['graph_summary_traversal'], path + '.traversal')
                check.ok(is_int(p['k']) and p['k'] >= 2, path + '.k')
                if p['available']:
                    check.ok(p['unavailable_reason'] is None, path + '.unavailable_reason')
                    check.ok(p['counting'] == ('upper_bound' if p['mask'] == 'absent'
                                               else 'exact'), path + '.counting')
                else:
                    check.one_of(p['unavailable_reason'], UNAVAILABLE, path + '.unavailable_reason')
                    check.ok(p['counting'] is None, path + '.counting')
                check.ok(p['mask'] is None or p['mask'] in MASKS, path + '.mask')
                check.ok(p['graph_mode'] is None or p['graph_mode'] in GRAPH_MODES,
                         path + '.graph_mode')
            check.ok(b['available'] is any(p['available'] for p in summary['pairs']),
                     'pattern.available', 'whether a pair is served')

    def test_the_block_is_the_same_on_both_routes(self):
        """SPEC §23: on every fixture server the full block is the same on /capabilities and
        /pattern/capabilities, and the block of /traverse/capabilities is the full block with
        `details` (phase 1), every gate field with the full block's value."""
        by_server = {}
        for name, f in self.fixtures.items():
            if f['method'] == 'GET' and f['status'] == 200 and not f['hand_made']:
                by_server.setdefault(f['server'], []).append(name)
        routes = set()
        for server, names in by_server.items():
            full = self.capabilities_of(server)
            for name in names:
                with self.subTest(fixture=name):
                    route = self.fixtures[name]['path'].split('?')[0]
                    routes.add(route)
                    if server == 'multi' and route == '/capabilities':
                        # a multi-graph server's own block: its pair's without a graph (SPEC
                        # §24), as pattern_capabilities_json writes it without one, and
                        # available since its pair is served
                        expected = copy.deepcopy(full)
                        for f in GRAPH_FIELDS:
                            expected[f] = None
                        expected['predicate']['access'] = None
                        self.assertEqual(expected, self.full_block(name))
                        continue
                    self.assertEqual(full, self.full_block(name))
                    if route != '/traverse/capabilities':
                        continue
                    b = self.bodies[name][1]['pattern']
                    self.assertEqual(dict(full, details=DETAILS), b)
                    for key in GATE_KEYS:
                        if key in full:
                            self.assertEqual(full[key], b[key], key)
                    for cap in GATE_CAPS if 'caps' in full else ():
                        self.assertEqual(full['caps'][cap], b['caps'][cap], cap)
        self.assertEqual({'/capabilities', '/traverse/capabilities', '/pattern/capabilities'},
                         routes)

    def test_gate_keys_are_the_servers(self):
        """GATE_KEYS and GATE_CAPS are the server's lists (pattern.cpp kPatternGateKeys,
        kPatternGateCaps, read from the source), each a field of the full block (SPEC §10.2)."""
        self.assertLessEqual(set(GATE_KEYS), set(SCHEMA['capabilities']))
        self.assertLessEqual(set(GATE_CAPS), set(CAPS))
        self.assertEqual(len(GATE_KEYS) + len(GATE_CAPS), len(set(GATE_KEYS) | set(GATE_CAPS)))
        if not os.path.isfile(PATTERN_CPP):
            self.skipTest('the server sources are not in this checkout')
        with open(PATTERN_CPP, encoding='utf-8') as f:
            source = f.read()
        for name, expected in (('kPatternGateKeys', GATE_KEYS), ('kPatternGateCaps', GATE_CAPS)):
            m = re.search(r'const char \*const ' + name + r'\[\] = \{(.*?)\};', source, re.S)
            self.assertTrue(m, name)
            self.assertEqual(expected, re.findall(r'"([a-z_]+)"', m.group(1)), name)
        self.assertIn('"' + DETAILS + '"', source)

    def test_caps_rule_section_classifies_every_cap(self):
        """caps_rule refers to SPEC §4.5, and §4.5 classifies every cap of every fixture
        server's block: each named once, as a request field's maximum or as the server's
        policy."""
        if not os.path.isfile(SPEC):
            self.skipTest('docs/SPEC-pattern-search.md is not in this checkout')
        with open(SPEC, encoding='utf-8') as f:
            text = f.read()
        m = re.search(r'^### 4\.5 .*?$(.*?)^## 5\. ', text, re.M | re.S)
        self.assertTrue(m, 'SPEC §4.5')
        section = m.group(1)
        m = re.search(r'\*\*Every cap of the capabilities\' `caps`, classified\*\*.*?\n\n'
                      r'(.*?)\n\n', section, re.S)
        self.assertTrue(m, 'the classification in §4.5')
        kinds = re.findall(r'^- \*\*(.*?)\*\*(.*?)(?=^- |\Z)', m.group(1), re.M | re.S)
        self.assertEqual(['the maximum of a request field', "the server's policy"],
                         [k for k, _ in kinds])
        named = [re.findall(r'`([a-z_]+)`', body) for _, body in kinds]
        caps = set()
        for name, f in self.fixtures.items():
            b = self.full_block(name) if f['method'] == 'GET' else None
            if b and 'caps' in b:
                self.assertEqual(CAPS_RULE, b['caps_rule'], name)
                caps |= set(b['caps'])
        self.assertEqual(set(CAPS), caps)
        for cap in caps:
            # in exactly one of the two lists (the request field's name is the cap's)
            self.assertEqual(1, sum(cap in n for n in named), cap)
        policy = {'max_patterns', 'min_information_bits', 'max_checked_entries',
                  'max_predicate_labels', 'row_cache_mb'}
        self.assertEqual(policy, set(named[1]) & caps)
        # a request field's maximum is a field of the request
        self.assertLessEqual(caps - policy, set(SCHEMA['request']))

    def test_capabilities_routes_refuse_what_v1_never_answers(self):
        """SPEC §23: the stored capabilities documents pass, and each rule refuses a document
        that breaks it -- the /traverse block without `details`, naming another route, or
        missing a gate field; `details` in the full block; /capabilities without the route of
        the full block; caps_rule not the reference to §4.5; a 400 of ?graph= with a code or on
        a route that answers it; and the full block while the index loads with a graph field
        set."""
        def run(name, mutate=None):
            f = self.fixtures[name]
            doc = copy.deepcopy(self.bodies[name][1])
            if mutate:
                mutate(doc)
            self.capabilities_document(Checker(Raising(), name), f, doc)

        def pop(*path):
            def m(doc):
                for key in path[:-1]:
                    doc = doc[key]
                doc.pop(path[-1])
            return m

        def put(value, *path):
            def m(doc):
                for key in path[:-1]:
                    doc = doc[key]
                doc[path[-1]] = value
            return m

        for name, f in self.fixtures.items():
            if f['method'] == 'GET':
                run(name)
        cases = [
            ('traverse_capabilities', pop('pattern', 'details'), 'details'),
            ('traverse_capabilities', put('GET /capabilities', 'pattern', 'details'), 'details'),
            ('traverse_capabilities_multi_graph', pop('pattern', 'details'), 'details'),
            ('traverse_capabilities', pop('pattern', 'counting'), 'gate field'),
            ('traverse_capabilities_mask_absent', pop('pattern', 'long_search'), 'gate field'),
            ('pattern_capabilities', put(DETAILS, 'details'), 'fields'),
            ('capabilities', put(DETAILS, 'pattern', 'details'), 'details'),
            ('capabilities', pop('routes', 'pattern_capabilities'), 'pattern_capabilities'),
            ('capabilities_multi_graph', put(DETAILS, 'routes', 'pattern_capabilities'),
             'pattern_capabilities'),
            ('pattern_capabilities', put('SPEC-pattern-search.md section 7.4', 'caps_rule'),
             'section 4.5'),
            ('pattern_capabilities_graph_param', put('invalid_request', 'code'), 'fields'),
            ('pattern_capabilities_loading', put('basic', 'graph_mode'), 'loads'),
            ('pattern_capabilities_loading', put(31, 'k'), 'loads'),
        ]
        for name, mutate, says in cases:
            with self.subTest(fixture=name, says=says):
                with self.assertRaisesRegex(AssertionError, says):
                    run(name, mutate)
        # a 400 where the route answers: ?graph= on the multi-graph server, none on a single one
        f = dict(self.fixtures['pattern_capabilities_graph_param'], server='multi')
        with self.assertRaisesRegex(AssertionError, 'only \\?graph='):
            self.capabilities_document(Checker(Raising(), 'multi'), f,
                                       self.bodies['pattern_capabilities_graph_param'][1])
        f = dict(self.fixtures['pattern_capabilities_multi_graph_no_graph'], server='masked')
        with self.assertRaisesRegex(AssertionError, 'only \\?graph='):
            self.capabilities_document(
                Checker(Raising(), 'single'), f,
                self.bodies['pattern_capabilities_multi_graph_no_graph'][1])

    def test_multi_graph_answers(self):
        """SPEC §24.1: a multi-graph server's answers are the pair answers, tagged, each with
        its outcome, a refused pair's entry its own refusal; each rule refuses an answer that
        breaks it -- a tag missing, another index_fp, the names unsorted or another count of
        answers, a field of the envelope more, an outcome unknown or miscounted, a refusal in
        an answered entry, an answer's field in a refused one, a refusal without its status,
        with an unknown code, or with the deadline's code under another status."""
        self.assertTrue(self.multi_answers, 'the multi-graph fixtures')

        def run(name, mutate=None):
            a = copy.deepcopy(self.multi_answers[name])
            if mutate:
                mutate(a)
            Checker(Raising(), name).multi_answer(a, self.bodies[name][0], 1)

        for name, a in self.multi_answers.items():
            run(name)
            # the answer the other checks read is the answered pair's without its tags
            answered = [x for x in a['answers'] if x['outcome'] == 'answered']
            pair = {k: v for k, v in answered[0].items() if k not in SCHEMA['answer_pair']}
            self.assertEqual(pair, self.bodies[name][1], name)
        name = 'multi_graph_count'
        cases = [
            (lambda a: a['answers'][0].pop('annotation_path'), 'tagged'),
            (lambda a: a['answers'][0].update(index_fp='0' * 64), 'index_fp'),
            (lambda a: a['graphs'].append('a_name'), 'byte order'),
            (lambda a: a['answers'].append(copy.deepcopy(a['answers'][0])), 'one per pair'),
            (lambda a: a.update(index={}), 'fields'),
            (lambda a: a['timing'].update(load_ms=0), 'fields'),
            (lambda a: a['answers'][0].pop('outcome'), 'outcome'),
            (lambda a: a['answers'][0].update(outcome='maybe'), 'outcome'),
            (lambda a: a.update(answered=0, refused=1), 'count the entries'),
            (lambda a: a['answers'][0].update(refusal={'http_status': 400, 'error': 'x'}),
             'no refusal'),
        ]
        for mutate, says in cases:
            with self.subTest(says=says):
                with self.assertRaisesRegex(AssertionError, says):
                    run(name, mutate)
        # the refused entry: the mixed server's fixture
        name = 'multi_graph_refused_pair'
        self.assertIn(name, self.multi_answers)
        refused = self.multi_answers[name]['answers'][0]
        self.assertEqual(('hashed', 'refused', 400, 'representation_unsupported'),
                         (refused['graph'], refused['outcome'],
                          refused['refusal']['http_status'], refused['refusal']['code']))
        self.assertEqual((1, 1), (self.multi_answers[name]['answered'],
                                  self.multi_answers[name]['refused']))
        cases = [
            (lambda a: a['answers'][0].update(patterns=[]), 'answers\\[0\\]'),
            (lambda a: a['answers'][0]['refusal'].pop('http_status'), 'http_status'),
            (lambda a: a['answers'][0]['refusal'].update(code='nope'), 'code'),
            (lambda a: a['answers'][0]['refusal'].update(http_status=503), '503 is the deadline'),
            (lambda a: a['answers'][0]['refusal'].update(code='deadline'), '503 is the deadline'),
            (lambda a: a['answers'][0]['refusal'].update(limits={}), 'http_status, error'),
            (lambda a: a.update(answered=2, refused=0), 'count the entries'),
        ]
        for mutate, says in cases:
            with self.subTest(fixture=name, says=says):
                with self.assertRaisesRegex(AssertionError, says):
                    run(name, mutate)

    def test_multi_graph_capabilities_refuse_what_v1_never_answers(self):
        """SPEC §24: GET /capabilities of a multi-graph server and its blocks per pair; each
        rule refuses a document that breaks it -- a graph field set in the server-wide block,
        availability that no pair has, a summary pair's counting against its mask, a pair
        missing from the summary, columns_disjoint against shared_columns, the routes without
        ?graph=, a pair's block not naming its pair, and the 400 without ?graph= carrying a
        code."""
        def run(name, mutate=None):
            doc = copy.deepcopy(self.bodies[name][1])
            if mutate:
                mutate(doc)
            self.capabilities_document(Checker(Raising(), name), self.fixtures[name], doc)

        def summary_pair(**fields):
            def m(doc):
                doc['graph_summary']['pairs'][0].update(fields)
            return m

        cases = [
            ('capabilities_multi_graph', lambda d: d['pattern'].update(k=31), 'per pair'),
            ('capabilities_multi_graph', summary_pair(available=False,
                                                      unavailable_reason='mask_invalid',
                                                      counting=None),
             'whether a pair is served'),
            ('capabilities_multi_graph', summary_pair(counting='upper_bound'), 'counting'),
            ('capabilities_multi_graph', lambda d: d['graph_summary']['pairs'].clear(),
             'every name'),
            ('capabilities_multi_graph', lambda d: d['graph_summary'].update(shared_columns=2),
             'disjoint exactly'),
            ('capabilities_multi_graph', summary_pair(traversal={}), 'fields'),
            ('capabilities_multi_graph',
             lambda d: d['routes'].update(pattern_capabilities=DETAILS), 'pattern_capabilities'),
            ('capabilities', lambda d: d.update(graph_summary={}), 'null on one graph'),
            ('pattern_capabilities_multi_graph', lambda d: d.pop('graph'), 'the pair named'),
            ('pattern_capabilities_multi_graph', lambda d: d.update(graph='x'), 'the name'),
            ('pattern_capabilities_multi_graph_no_graph',
             lambda d: d.update(code='invalid_request'), 'fields'),
        ]
        for name, mutate, says in cases:
            with self.subTest(fixture=name, says=says):
                run(name)
                with self.assertRaisesRegex(AssertionError, says):
                    run(name, mutate)

    def test_requests_and_answers_agree(self):
        """What a fixture's request asks is what its answer states."""
        for name, f in self.fixtures.items():
            request, answer = self.bodies[name]
            if f['method'] != 'POST' or f['status'] != 200:
                continue
            with self.subTest(fixture=name):
                # without the mask the thresholds compare the upper bound (conservative), and
                # the note threshold_upper_bound says when the lower bound was within the
                # threshold
                unmasked = 'counting' in answer['index']
                threshold = answer['limits']['max_contexts']
                for e in answer['patterns']:
                    if 'error' in e or e['scope'] == 'long':
                        continue
                    noted = 'threshold_upper_bound' in e['notes']
                    if request.get('stop_at_threshold') and e['stop'] \
                            and e['stop']['reason'] == 'max_contexts' \
                            and e['stop']['phase'] != 'selection':
                        c = e['counts']['contexts']
                        self.assertEqual('at_least', c['relation'])
                        if unmasked:
                            self.assertTrue(c['value'] > threshold or noted, c)
                        else:
                            self.assertGreater(c['value'], threshold)
                    if e.get('withheld') == {'reason': 'count_above_threshold'}:
                        c = e['counts']['contexts']
                        if unmasked and c['relation'] == 'bounds':
                            self.assertGreater(c['upper'], threshold)
                            self.assertEqual(c['lower'] <= threshold, noted, c)
                        else:
                            self.assertEqual('exact', c['relation'])
                            self.assertGreater(c['value'], threshold)
                    if not unmasked:
                        self.assertFalse(noted, 'threshold_upper_bound with the mask')
                    if e.get('cut') == {'reason': 'max_contexts'}:
                        self.assertEqual(threshold, e['returned'])

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
        # paths and peptides: every kind; the paths' extension completed, stopped, not
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
        self.assertEqual({'max_contexts', 'max_steps', 'time', 'max_memory', 'max_paths',
                          'max_predicate_contexts', 'max_predicate_work',
                          'max_annotation_work'}, seen['cut'])
        # a predicate: the selection's stops and the raw threshold's
        self.assertLessEqual({('selection', 'max_predicate_work'), ('selection', 'max_contexts'),
                              ('discovery', 'max_predicate_contexts')}, seen['stop'])
        self.assertEqual(set(SLOT_ERRORS), seen['slot'])
        self.assertEqual({'basic', 'primary'}, seen['graph_mode'])
        self.assertEqual(set(NOTES), seen['note'])
        # every relation, bounds and the mask_scan stop included
        self.assertEqual(set(RELATIONS), seen['relation'])
        self.assertLessEqual({('discovery', 'max_steps'), ('discovery', 'time'),
                              ('discovery', 'max_contexts'), ('mask_scan', 'max_steps')},
                             seen['stop'])
        # every refusal code but the unproducible ones (a set pinned to the covered ones would
        # let a code without a fixture go unseen) and the ones no stored fixture holds, listed
        # by name with the index each needs (derived from REFUSALS alone, the expectation could
        # not notice REFUSALS missing codes; test_the_codes_are_the_sources compares it with the
        # code); the retired mask_required has none either (RETIRED)
        without_fixture = {'mask_invalid', 'alphabet_untested'}
        self.assertEqual(without_fixture, set(NO_FIXTURE))
        self.assertEqual({'mask_required', 'resident_only'}, set(RETIRED))
        unproducible = set(UNPRODUCIBLE) & set(REFUSALS)
        self.assertEqual((set(REFUSALS) | {'initializing'}) - unproducible
                         - without_fixture - set(RETIRED), seen['refusal'])
        self.assertEqual(without_fixture | unproducible | set(RETIRED),
                         set(REFUSALS) - seen['refusal'], 'the refusal codes without a fixture')
        self.assertLessEqual({('label_discovery', 'max_annotation_work')}, seen['stop'])
        # the extension's stops
        self.assertLessEqual({('extension', 'max_paths'), ('extension', 'max_steps')},
                             seen['stop'])
        # SPEC §7.6: stop keeps the first stop; the time that then cut the release shows only
        # as cut time and time_limited, in the pattern stopped by max_steps and in the one after
        # it
        self.assertIn(('max_steps', 'time', 'time_limited'), seen['stop_then_cut'])
        # supported paths (SPEC §20): what the search did, its stops (the work budget, a row
        # the account refused, the list's memory), both levels, a predicate's pruning and the
        # mirror walks read; the motif (SPEC §25): each way it is decided, and undecided
        searches, stops, levels, phases, motifs = set(), set(), set(), set(), set()
        pruned = mirrors = 0
        for name, f in self.fixtures.items():
            request, answer = self.bodies[name]
            if f['method'] != 'POST' or f['status'] != 200:
                continue
            for e in answer['patterns']:
                if 'supported_paths' in e.get('counts', {}):
                    c = e['counts']['supported_paths']
                    searches.add(c['search'])
                    levels.add(c['level'])
                    pruned += c.get('branches_pruned_by_predicate', 0) > 0
                    mirrors += e['work'].get('mirror_rows', 0) > 0
                    if e['stop']:
                        stops.add((e['stop']['phase'], e['stop']['reason']))
                    phases.update(x['phase'] for x in e.get('rows_refused') or [])
                if 'motif' in e:
                    motifs.add((e['motif']['decided_by'], e['motif']['untested']))
        self.assertLessEqual({'completed', 'stopped', 'no_anchors'}, searches)
        self.assertEqual({'record_verified', 'label_intersection'}, levels)
        self.assertLessEqual({('extension', 'max_annotation_work'), ('extension', 'max_memory'),
                              ('output', 'max_memory')}, stops)
        self.assertIn('extension', phases)
        self.assertTrue(pruned and mirrors)
        self.assertLessEqual({('every_context', None), ('tested_contexts', 'selection'),
                              ('no_instance', None), (None, 'not_started')}, motifs)

    def test_time_limited_is_the_clocks(self):
        """SPEC §7.6/§7.9, both ways: an entry the clock touched (its stop, or only its cut) is
        time_limited, and a time_limited entry was touched by the clock (the SPEC's stop
        max_steps + cut time is accepted)."""
        name = 'max_steps_then_time'
        request, answer = self.bodies[name]
        capabilities = self.capabilities_of(self.fixtures[name]['server'])

        def check(a):
            Checker(Raising(), name).answer(a, request, capabilities)

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

    def test_every_capability_value_has_a_capabilities_fixture(self):
        """On both capabilities routes: every unavailable reason but the unproducible ones and
        the ones named as without a fixture (NO_FIXTURE); every counting (exact with a mask file
        or built at load, upper_bound without, with its dummy_fraction; the answers of each
        server state the same, checked in Checker.answer); and every mask value (file,
        built_at_load and absent, DESIGN §4)."""
        seen = {'unavailable': {}, 'counting': {}, 'mask': {}}
        for name, f in self.fixtures.items():
            if f['method'] != 'GET' or f['status'] != 200 or f['hand_made']:
                continue
            route = 'probe' if f['path'].startswith('/traverse/') \
                else 'pattern' if f['path'].startswith('/pattern/') else 'capabilities'
            b = self.full_block(name)
            if b['available'] is False:
                seen['unavailable'].setdefault(b['unavailable_reason'], set()).add(route)
            if b.get('counting') is not None:
                seen['counting'].setdefault(b['counting'], set()).add(route)
            # (a multi-graph server and a graph the engine does not recognise state no mask)
            if b.get('mask') is not None:
                seen['mask'].setdefault(b['mask'], set()).add(route)
        without_fixture = {'mask_invalid', 'alphabet_untested'}
        self.assertEqual(without_fixture, set(NO_FIXTURE))
        self.assertEqual(set(UNAVAILABLE) - set(UNPRODUCIBLE) - without_fixture
                         - RETIRED_REASONS,
                         set(seen['unavailable']))
        self.assertEqual(without_fixture | set(UNPRODUCIBLE) | RETIRED_REASONS,
                         set(UNAVAILABLE) - set(seen['unavailable']),
                         'the unavailable reasons without a fixture')
        self.assertEqual(set(COUNTINGS), set(seen['counting']))
        self.assertEqual(set(MASKS), set(seen['mask']))
        for kind, values in seen.items():
            for value, routes in values.items():
                self.assertLessEqual({'probe', 'capabilities'}, routes, (kind, value))

    def test_the_codes_are_the_sources(self):
        """REFUSALS and UNAVAILABLE are the codes the server's sources write (a code missing
        here would go unnoticed by the coverage tests, which are derived from these lists):
        every PatternRefusal's literal code, and every reason support_message words (the
        graph-support reasons a refusal and the capabilities' unavailable_reason carry), with
        the multi-graph server's own unavailable reason."""
        if not os.path.isfile(PATTERN_CPP):
            self.skipTest('the server sources are not in this checkout')
        cli = os.path.dirname(PATTERN_CPP)
        sources = {}
        for name in sorted(os.listdir(cli)):
            if name.endswith(('.cpp', '.hpp')):
                with open(os.path.join(cli, name), encoding='utf-8') as f:
                    sources[name] = f.read()
        literal = set()
        for name, text in sources.items():
            codes = set(re.findall(r'PatternRefusal\(\s*\d+\s*,\s*"(\w+)"', text))
            if name in NOT_SERVED_SOURCES:
                # a module that is built and tested but not called by the route (the
                # request field it reads is refused by name, later_increment): its codes are
                # no answer of this build
                self.assertEqual(set(NOT_SERVED_SOURCES[name]), codes - set(REFUSALS), name)
                self.assertNotIn(f'#include "{name[:-len(".cpp")]}.hpp"', sources['pattern.cpp'],
                                 f'{name} is served now: its codes belong in REFUSALS')
                continue
            literal |= codes
        support = set(re.findall(r'support\.reason == "(\w+)"', sources['pattern.cpp']))
        self.assertIn('PatternRefusal(400, support.reason, support_message(support))',
                      sources['pattern.cpp'], 'a graph-support refusal is its reason')
        unavailable = set(re.findall(r'p\["unavailable_reason"\] = "(\w+)"',
                                     sources['pattern.cpp']))
        self.assertLessEqual({'mask_invalid', 'alphabet_untested'}, support,
                             'support_message\'s branches not found in pattern.cpp')
        # the retired codes (mask_required) are written by no source: a code is written as a
        # support reason (reason = "..." in the engine and the route) or a literal refusal
        engine = os.path.join(REPO, 'src', 'graph', 'alignment', 'pattern_search.cpp')
        with open(engine, encoding='utf-8') as f:
            engine_text = f.read()
        assigned = set(re.findall(r'reason = "(\w+)"', engine_text + sources['pattern.cpp']))
        self.assertFalse(set(RETIRED) & (literal | support | unavailable | assigned),
                         'a retired code is still written')
        self.assertEqual(set(REFUSALS), literal | support | set(RETIRED))
        self.assertLessEqual(RETIRED_REASONS, set(RETIRED))
        self.assertEqual(set(UNAVAILABLE), support | unavailable | RETIRED_REASONS)
        # the slot codes: the engine's PatternError code (bad_alphabet; stop_unsupported is
        # retired) and its refusals of a parsed pattern (information_below_floor,
        # scope_unsupported)
        slots = set(re.findall(r'PatternError\(\s*"(\w+)"', engine_text))
        self.assertEqual({'bad_alphabet'}, slots)
        self.assertNotIn('stop_unsupported', engine_text)
        for code in ('information_below_floor', 'scope_unsupported'):
            self.assertIn(f'"{code}"', engine_text)
        self.assertEqual(set(SLOT_ERRORS), slots | {'information_below_floor',
                                                    'scope_unsupported'})

    def test_by_label_null_in_partial_is_valid(self):
        """SPEC §14.4: in partial, when the memory account cannot hold by_label, the answer has
        by_label null, every context read output_budget and a stop -- a valid v1 answer the
        validator must accept. The body is a real CLI answer on a tiny long-label index
        (data/traverse/pattern_validator/by_label_null_partial); no fixture of the mini index
        can show it (its label names are short)."""
        d = os.path.join(HERE, 'data', 'traverse', 'pattern_validator', 'by_label_null_partial')
        with open(os.path.join(d, 'request.json')) as f:
            request = json.load(f)
        # stored gzipped: one label name of 512 KiB
        with gzip.open(os.path.join(d, 'answer.json.gz'), 'rt') as f:
            answer = json.load(f)
        e = answer['patterns'][1]
        self.assertEqual((None, {'phase': 'output', 'reason': 'max_memory'}, 'exact'),
                         (e['by_label'], e['stop'], e['counts']['labels']['relation']))

        def check(a, req=request):
            Checker(Raising(), 'by_label_null_partial').answer(a, req)

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
        accepted (none is stored, see NO_FIXTURE), and the rules they rest on still refuse what
        v1 never answers."""
        def refusal(code, status=400):
            check = Checker(Raising(), f'hand-made {code}')
            self.refusal(check, {'status': status, 'headers': {}},
                         {'error': f'pattern: the graph is not served ({code})', 'code': code})

        def block(**fields):
            b = copy.deepcopy(self.capabilities_of('masked'))
            self.assertIsNotNone(b)
            self.assertIs(True, b['available'])
            b.update(fields)
            if b['available'] is not True and 'counting' not in fields:
                # a graph not served states no counting
                b['counting'] = None
            Checker(Raising(), 'hand-made capabilities').block(b, 'pattern', False)

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
        with self.assertRaisesRegex(AssertionError, 'mask_invalid: a mask file'):
            block(available=False, unavailable_reason='mask_invalid', mask='absent')
        with self.assertRaisesRegex(AssertionError, 'is served'):
            block(available=False, unavailable_reason='alphabet_untested')
        with self.assertRaisesRegex(AssertionError, 'only \\$ACGT is served'):
            block(alphabet='$ACGTN')
        with self.assertRaisesRegex(AssertionError, 'unavailable_reason'):
            block(available=False, unavailable_reason='alphabet_untested', alphabet='$ACGU')
        # mask_required is retired (an older build's block, which has no counting, is not this
        # build's); counting and the dummy fraction follow the mask
        with self.assertRaisesRegex(AssertionError, 'a retired reason'):
            block(available=False, unavailable_reason='mask_required', mask='absent')
        with self.assertRaisesRegex(AssertionError, 'exact with a mask, upper_bound without'):
            block(mask='absent')
        with self.assertRaisesRegex(AssertionError, 'exact with a mask, upper_bound without'):
            block(counting='upper_bound')
        with self.assertRaisesRegex(AssertionError, 'null on a graph not served'):
            block(available=False, unavailable_reason='mask_invalid', counting='exact')
        unmasked = self.capabilities_of('unmasked')
        with self.assertRaisesRegex(AssertionError, 'null with the mask'):
            block(dummy_fraction=unmasked['dummy_fraction'])
        with self.assertRaisesRegex(AssertionError, 'an interval in \\[0, 1\\] holding'):
            block(mask='absent', counting='upper_bound',
                  dummy_fraction=dict(unmasked['dummy_fraction'], interval=[0.5, 0.6]))
        # a v1 refusal mask_required, as an older build of version 1 answers it, is still a
        # refusal a client handles (the code stays in the contract's vocabulary)
        refusal('mask_required')

    def test_paths_and_peptides_rules_refuse_what_v1_never_answers(self):
        """Paths and peptides (SPEC §12.1, §12.2): the stored path and peptide answers pass,
        and each rule they rest on refuses an answer that breaks it (so that the rules are not
        vacuous): a path with kmer, a path that does not instantiate the pattern, a label
        stated record_verified without an occurrence, an unverified label listed under
        require_support, a stopped extension stated exact, a peptide's bits, a peptide's
        instance off its codons, a stop refused as bad_alphabet, a path's excluded count stated
        over truncated rows or withheld over complete ones."""
        for name in ('paths', 'paths_labels', 'paths_require_support', 'paths_stop_at_max_paths',
                     'peptide', 'peptide_paths', 'peptide_bad_residue', 'paths_global',
                     'paths_primary'):
            self.check_stored(name)

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
            # the stop '*' is a residue: refusing it is not v1's answer
            ('peptide_bad_residue',
             lambda a: a['patterns'].__setitem__(1, {
                 'id': 'stop', 'kind': 'protein',
                 'error': {'code': 'bad_alphabet', 'message': 'pattern: *'}}),
             'bad_alphabet for a pattern inside the alphabet'),
            # a path's labels_excluded_unverified: an integer for a path whose rows were
            # truncated, and null for a path read and verified completely in an answer whose
            # count is exact
            ('paths_require_support',
             lambda a: first_path(a).update(labels_status='truncated', labels_total=None),
             'null unless every row of the path was read completely'),
            ('paths_require_support',
             lambda a: a['patterns'][1]['results'][0].update(labels_excluded_unverified=None),
             'undecided only when the verification was not done'),
        ]
        self.assertRefusesMutations(cases)

    def test_hand_made_bodies_are_the_codes(self):
        """The hand-made bodies are written by the code as stored here: the two 503s, and the
        full block while the index loads (the loaded masked server's block with what
        pattern_capabilities_json sets to null without the graph)."""
        if not (os.path.isfile(PATTERN_CPP) and os.path.isfile(SERVER_UTILS_CPP)):
            self.skipTest('the server sources are not in this checkout')
        with open(PATTERN_CPP, encoding='utf-8') as f:
            pattern_cpp = f.read()
        with open(SERVER_UTILS_CPP, encoding='utf-8') as f:
            server_utils = f.read()
        hand_made = {n for n, f in self.fixtures.items() if f['hand_made']}
        self.assertEqual({'deadline_503', 'initializing_503', 'pattern_capabilities_loading'},
                         hand_made)

        # the graph fields the code sets to null without the graph, read from its source
        m = re.search(r'const char \*graph_fields\[\] = \{(.*?)\};', pattern_cpp, re.S)
        self.assertTrue(m, 'pattern_capabilities_json graph_fields')
        self.assertEqual(list(GRAPH_FIELDS), re.findall(r'"([a-z_]+)"', m.group(1)))
        loading = self.bodies['pattern_capabilities_loading'][1]
        loaded = self.bodies['pattern_capabilities'][1]
        self.assertTrue(loaded['available'])
        self.assertEqual(set(loaded), set(loading))
        for key, value in loading.items():
            if key in GRAPH_FIELDS + ('available', 'unavailable_reason'):
                self.assertIsNone(value, key)
            elif key == 'predicate':
                self.assertEqual(dict(loaded['predicate'], access=None), value)
            else:
                self.assertEqual(loaded[key], value, key)

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

    def test_mask_required_is_retired(self):
        """A graph without its dummy-edge mask is served, so no fixture holds mask_required (the
        unmasked_* fixtures show such graphs), support_message has no branch for it, and the
        unmasked server's capabilities say available, counting upper_bound (the codes against
        the sources: test_the_codes_are_the_sources)."""
        self.assertNotIn('mask_required', self.fixtures)
        for name, (request, answer) in self.bodies.items():
            self.assertNotIn('"mask_required"', json.dumps(answer), name)
        b = self.capabilities_of('unmasked')
        self.assertEqual((True, None, 'absent', 'upper_bound'),
                         (b['available'], b['unavailable_reason'], b['mask'], b['counting']))
        unmasked = [n for n, f in self.fixtures.items()
                    if f['server'] == 'unmasked' and f['method'] == 'POST']
        self.assertTrue(unmasked and all(self.fixtures[n]['status'] == 200 for n in unmasked))
        if os.path.isfile(PATTERN_CPP):
            with open(PATTERN_CPP, encoding='utf-8') as f:
                self.assertNotIn('support.reason == "mask_required"', f.read())

    def test_unmasked_and_stop_rules_refuse_what_v1_never_answers(self):
        """Graphs without their mask and the stop '*': the stored answers without the mask and
        with the stop '*' pass, and each rule they rest on refuses an answer that breaks it: an
        estimate off its formula, a bounds count without its estimate, an estimate on a masked
        graph, the estimate note missing, a bounds count stated without a stop on a masked
        graph, threshold_upper_bound missing or on a masked graph, an index without its dummy
        fraction, a '*' instance that is no stop codon, no_stop_codon missing or a no-instance
        peptide with a context, and the unsearched peptide after a stop stated stopped; and
        bounds left on a pattern with no more unchecked candidates than the server's
        max_checked_entries."""
        def total(a, i=0):
            c = a['patterns'][i]['counts']
            return c.get('contexts') or c['anchors']

        for name in ('unmasked_count', 'unmasked_labels_all', 'unmasked_threshold_upper_bound',
                     'unmasked_stop_at_threshold', 'unmasked_partial', 'unmasked_paths',
                     'unmasked_checked', 'unmasked_checked_dummies',
                     'peptide_stop', 'peptide_no_stop_codon', 'peptide_no_stop_codon_after_stop',
                     'peptide_bad_residue'):
            self.check_stored(name)

        def drop_note(i, note):
            return lambda a: a['patterns'][i]['notes'].remove(note)

        def masked(a):
            # the same body as if it came from a masked graph
            for f in ('counting', 'dummy_fraction'):
                del a['index'][f]

        def unchecked_start(a):
            # the island start's entry of unmasked_count (the check off), in unmasked_checked
            e = self.bodies['unmasked_count'][1]['patterns'][1]
            self.assertEqual('start', e['id'])
            a['patterns'][1].update({f: copy.deepcopy(e[f]) for f in ('counts', 'notes', 'work')})

        def stop_codon_read_as(codon):
            # the first context's instance (and its k-mer) with its stop codon TGA replaced
            def run(a):
                r = a['patterns'][0]['results'][0]
                instance = r['instance'][:-3] + codon
                o = r['offset']
                r.update(instance=instance,
                         kmer=r['kmer'][:o] + instance + r['kmer'][o + len(instance):])
            return run

        cases = [
            ('unmasked_count', lambda a: total(a, 1).update(estimate=17), 'estimate'),
            ('unmasked_count', lambda a: total(a).pop('estimate'), 'fields'),
            ('unmasked_count', lambda a: total(a, 2).update(estimate=0), 'fields'),
            ('unmasked_count', drop_note(0, 'estimate_sampled_dummy_fraction'),
             'estimate_sampled_dummy_fraction exactly'),
            ('unmasked_count', masked, 'fields|counting'),
            ('unmasked_count', lambda a: a['index'].pop('dummy_fraction'), 'fields'),
            ('unmasked_count', lambda a: a['index']['dummy_fraction'].update(source='guessed'),
             'source'),
            ('unmasked_threshold_upper_bound', drop_note(0, 'threshold_upper_bound'),
             'threshold_upper_bound exactly where the lower bound'),
            ('unmasked_stop_at_threshold', drop_note(0, 'threshold_upper_bound'),
             'a threshold stop below the threshold'),
            # TGG is W, not a stop: '*' admits the table's stop codons only
            ('peptide_stop', stop_codon_read_as('TGG'), 'does not instantiate'),
            ('peptide_no_stop_codon', drop_note(0, 'no_stop_codon'), 'no_stop_codon exactly'),
            ('peptide_no_stop_codon', lambda a: total(a).update(value=1),
             'no instance: exact 0|sum'),
            ('peptide_no_stop_codon_after_stop',
             lambda a: a['patterns'][1]['work'].update(mask_scans=1),
             'answered without a search'),
            # the island start answered as with the check off (its 30 unchecked candidates in
            # bounds [2, 32]) by the server that checks 50
            ('unmasked_checked', unchecked_start, 'they are checked, exact'),
        ]
        self.assertRefusesMutations(cases)
        # a masked graph never states bounds without a stop, nor threshold_upper_bound
        request, answer = self.bodies['count']
        a = copy.deepcopy(answer)
        c = a['patterns'][0]['counts']['contexts']
        for x in [c, c['suffix']] + list(c['by_offset'].values()) \
                + list(c['by_strand'].values()):
            x.update(relation='bounds', lower=x['value'], upper=x['value'])
        with self.assertRaisesRegex(AssertionError, 'without a stop is exact'):
            Checker(Raising(), 'count').answer(a, request, self.capabilities_of('masked'))
        a = copy.deepcopy(answer)
        a['patterns'][0]['notes'].append('threshold_upper_bound')
        with self.assertRaisesRegex(AssertionError, 'only on a graph without its mask'):
            Checker(Raising(), 'count').answer(a, request, self.capabilities_of('masked'))

    def test_round_fix3_rules_refuse_what_v1_never_answers(self):
        """The work counters, the low-complexity note and the capabilities' references: the
        stored bodies pass, and each rule refuses an answer that breaks it: a counter missing
        where its entry has it, or present where it has not; a completed extension that did not
        begin at every anchor; an extension that did not run with work of its own; more
        branchings than half the candidates entered; more distinct rows than rows read;
        verification steps without coordinates; low_complexity_pattern beside a stop;
        time_limited without a stop on a pattern its diagnostic reads whole; and capabilities
        whose prose fields are not references in ASCII, or whose delivery rates are not
        numbers."""
        def work(i, **fields):
            return lambda a: a['patterns'][i]['work'].update(fields)

        def drop(i, block, field):
            return lambda a: a['patterns'][i][block].pop(field)

        def with_placement(placement):
            def run(a):
                a['patterns'][0]['placement'] = placement
            return run

        for name in ('paths', 'paths_labels', 'paths_global', 'paths_anchors_above_threshold',
                     'labels_all', 'count_low_complexity', 'count'):
            self.check_stored(name)
        cases = [
            ('paths', drop(0, 'work', 'extension_anchors'), 'fields'),
            ('paths', drop(0, 'work', 'extension_branches'), 'fields'),
            ('count', work(0, extension_anchors=0), 'fields'),
            ('paths', work(0, extension_anchors=1), 'began at every anchor'),
            ('paths', work(0, extension_branches=10), 'two candidates per branching'),
            ('paths_anchors_above_threshold', work(0, extension_anchors=1), 'did not run'),
            ('paths', work(1, extension_branches=1), 'did not run'),
            ('labels_all', drop(0, 'work', 'annotation_rows_distinct'), 'fields'),
            ('labels_all', work(0, annotation_rows_distinct=49), 'a row read'),
            ('labels_all', work(0, verification_steps=0), 'fields'),
            ('paths_labels', drop(0, 'work', 'verification_steps'), 'fields'),
            ('paths_labels', drop(0, 'timing', 'verification_ms'), 'fields'),
            ('paths_labels', drop(0, 'timing', 'label_intersection_ms'), 'fields'),
            ('labels_all', lambda a: a['patterns'][0]['timing'].update(verification_ms=0.1),
             'fields'),
            ('paths_global',
             lambda a: (with_placement('none')(a),
                        a['patterns'][0]['notes'].__setitem__(0, 'label_intersection_only')),
             'nothing verified'),
            ('count_low_complexity',
             lambda a: a['patterns'][0].update(stop={'phase': 'extraction', 'reason': 'time'},
                                               determinism='time_limited'),
             'low_complexity_pattern beside a budget stop'),
            ('count', lambda a: a['patterns'][0].update(determinism='time_limited'),
             'time_limited without a time stop'),
        ]
        self.assertRefusesMutations(cases)
        # the capabilities
        for field, value, says in (('caps_rule', 'max_contexts ... see the SPEC', 'reference'),
                                   ('protein_rule', 'SPEC-pattern-search.md section §12.2',
                                    'ASCII'),
                                   ('delivery_mbps', {'build': 10}, 'fields'),
                                   ('delivery_mbps', {'build': 10, 'compress': '50'}, 'MB/s'),
                                   ('delivery_mbps', {'build': 0, 'compress': 50}, 'MB/s')):
            b = copy.deepcopy(self.capabilities_of('masked'))
            b[field] = value
            with self.subTest(field=field, value=value):
                with self.assertRaisesRegex(AssertionError, says):
                    Checker(Raising(), 'capabilities').block(b, 'pattern', False)

    def test_predicate_rules_refuse_what_v1_never_answers(self):
        """SPEC §19: the stored bodies pass, and each rule refuses an answer that breaks it --
        the predicate block missing or where no predicate was asked, a normal form that is not
        the request folded with its unknown names, a wrong vacuous, strands other than
        evaluated; a selection pass whose counts break §19.7's relations (completed but not
        every raw context tested, stopped without the bounds [S, S + R - T], constant tested
        other than the raw count, a selected count above the tested one); a complete answer
        whose results are not the selected count; withheld and cut reasons without their cause;
        a result's selection_labels the normal form does not hold on, or out of label order,
        its selection_strands missing, of another length, or naming a row its own labels
        contradict; predicate_only labels outside the predicate; the notes predicate_constant
        and projection_not_read where they do not apply, annotation_not_read on a predicate
        answer; the selection's fields in an answer without a predicate; and capabilities whose
        predicate object is not this build's."""
        def entry(i, **fields):
            return lambda a: a['patterns'][i].update(fields)

        def counts(i, name, **fields):
            return lambda a: a['patterns'][i]['counts'][name].update(fields)

        def block(**fields):
            return lambda a: a['predicate'].update(fields)

        def results(i, update):
            def run(a):
                for r in a['patterns'][i]['results']:
                    update(r)
            return run

        names = [n for n in self.fixtures if n.startswith('predicate_')
                 and self.fixtures[n]['status'] == 200]
        self.assertGreaterEqual(len(names), 19)
        for name in names:
            self.check_stored(name)
        cases = [
            ('predicate_filter', lambda a: a.pop('predicate'), 'fields'),
            ('count', lambda a: a.update(predicate=None), 'fields'),
            ('predicate_unknown_folded', block(normal_form={'any': ['562', '5622']}),
             'the request folded'),
            ('predicate_unknown_folded', block(unknown_labels=[]), 'names - unknown'),
            ('predicate_either', block(vacuous=False), 'vacuous'),
            ('predicate_context', block(strands='either'), 'evaluated'),
            ('predicate_unbudgeted_allowed', block(strands='context'), 'evaluated'),
            ('predicate_filter', counts(0, 'tested', value=1827), 'every raw context tested'),
            ('predicate_budget', counts(0, 'selected', upper=1828), 'bounds'),
            ('predicate_budget', counts(0, 'selected', relation='at_least'), 'fields'),
            ('predicate_unknown_constant', counts(0, 'tested', value=23), 'the raw count'),
            ('predicate_unknown_constant', counts(0, 'selected', value=24),
             'selected <= tested|0, or the raw count'),
            ('predicate_above_threshold', counts(0, 'selected', value=3, relation='exact'),
             'not_admitted: unknown'),
            ('predicate_filter', entry(0, returned=59,
                                       results=[]), 'the length of results|every selected'),
            ('predicate_filter', lambda a: a['patterns'][0]['results'].pop(),
             'every selected context returned|the length'),
            ('predicate_selected_above', counts(0, 'selected', value=5),
             'above max_contexts'),
            ('predicate_above_threshold',
             lambda a: a['patterns'][0]['selection'].update({'pass': 'not_started'}),
             'not admitted|not_started'),
            ('predicate_budget_partial', entry(0, cut={'reason': 'max_predicate_contexts'}),
             'the raw release held'),
            ('predicate_filter', results(0, lambda r: r.update(selection_labels=[])),
             'satisfy the normal form'),
            ('predicate_filter', results(0, lambda r: r.update(selection_labels=['562', '287'])),
             'satisfy the normal form|names of the normal form'),
            ('predicate_all', results(0, lambda r: r.update(selection_labels=['546'])),
             'satisfy the normal form'),
            ('predicate_only_record', results(0, lambda r: r['labels'][0].update(column='573')),
             'predicate_only|labels of by_label'),
            ('predicate_filter', results(0, lambda r: r.pop('selection_labels')), 'fields'),
            # per selected result and label its orientation (selection_strands)
            ('predicate_filter', results(0, lambda r: r.pop('selection_strands')), 'fields'),
            ('predicate_selection_strands',
             results(0, lambda r: r.update(selection_strands=['context'] * len(
                 r['selection_labels']))), 'on the own row exactly'),
            ('predicate_selection_strands',
             results(0, lambda r: r.update(selection_strands=['either'] * len(
                 r['selection_labels']))), 'states the row'),
            ('predicate_only_record', results(0, lambda r: r.update(selection_strands=[])),
             'one per selection label'),

            ('predicate_unknown_folded', lambda a: a['patterns'][0]['notes'].append(
                'predicate_constant'), 'predicate_constant'),
            ('predicate_either', lambda a: a['patterns'][0].update(notes=[]),
             'projection_not_read'),
            ('predicate_context', lambda a: a['patterns'][0]['notes'].append(
                'annotation_not_read'), 'annotation_not_read'),
            ('count', lambda a: a['patterns'][0].update(absence_filter='predicate'), 'fields'),
            ('count', lambda a: a['patterns'][0]['counts'].update(
                tested={'value': 24, 'relation': 'exact', 'unit': 'graph_contexts'}), 'fields'),
            ('predicate_long_anchors', lambda a: a['patterns'][0]['selection'].update(
                {'pass': 'completed'}), 'not_started'),
            ('predicate_stop_at_threshold', entry(0, stop={'phase': 'selection',
                                                           'reason': 'max_predicate_contexts'}),
             'max_predicate_contexts stops the raw discovery'),
            ('predicate_filter', lambda a: a['patterns'][0]['work'].update(predicate_lookups=0)
             or a['patterns'][0]['work'].pop('predicate_rows'), 'fields'),
            ('predicate_unbudgeted_allowed',
             lambda a: a['patterns'][0]['notes'].remove('annotation_unbudgeted'),
             'annotation_unbudgeted exactly'),
            ('predicate_filter', lambda a: a['limits'].update(predicate_strands='context'),
             'as requested'),
            ('predicate_filter', lambda a: a['limits'].pop('max_predicate_labels'), 'fields'),
        ]
        self.assertRefusesMutations(cases)
        # the capabilities (§19.12)
        for field, value, says in (('operators', ['any', 'all'], 'operators'),
                                   ('strands', ['either'], 'strands'),
                                   ('access', 'cells', 'access')):
            b = copy.deepcopy(self.capabilities_of('masked'))
            b['predicate'][field] = value
            with self.subTest(field=field, value=value):
                with self.assertRaisesRegex(AssertionError, says):
                    Checker(Raising(), 'capabilities').block(b, 'pattern', False)
        b = copy.deepcopy(self.capabilities_of('masked'))
        b['projections_later_increment'] = ['predicate_only']
        with self.assertRaisesRegex(AssertionError, 'projections'):
            Checker(Raising(), 'capabilities').block(b, 'pattern', False)

    def test_supported_paths_rules_refuse_what_v1_never_answers(self):
        """SPEC §20: the stored bodies pass, and each rule refuses an answer that breaks it --
        counts.paths not a plain count, or below the supported paths, or exact although a
        predicate pruned walks; counts.supported_paths with a relation its search does not give,
        a level other than the one asked for or the index's best, a predicate's pruning under
        "either"; work counters that are not rows read, mirror walks read where none is, a
        selection that reads rows of its own; a stop in the search that leaves the supported
        paths exact; the release reasons without their cause; selection_strands contradicting a
        path's own support; the notes paths_later_increment, projection_not_read,
        label_intersection_only and low_complexity_pattern where they do not apply; the timing,
        the limits and the placement of a supported-path answer; and requests the route refuses."""
        def entry(i, **fields):
            return lambda a: a['patterns'][i].update(fields)

        def counts(i, name, **fields):
            return lambda a: a['patterns'][i]['counts'][name].update(fields)

        def work(i, **fields):
            return lambda a: a['patterns'][i]['work'].update(fields)

        def results(i, update):
            def run(a):
                for r in a['patterns'][i]['results']:
                    update(r)
            return run

        def both(*mutations):
            def run(a):
                for m in mutations:
                    m(a)
            return run

        def at_least(i):
            # the supported paths and each orientation's at_least (their sums kept)
            def run(a):
                c = a['patterns'][i]['counts']['supported_paths']
                c['relation'] = 'at_least'
                for part in c['by_strand'].values():
                    part['relation'] = 'at_least'
            return run

        names = [n for n in self.fixtures if n.startswith('supported_paths')
                 and self.fixtures[n]['status'] == 200]
        self.assertGreaterEqual(len(names), 18)
        for name in names:
            self.check_stored(name)
        cases = [
            ('supported_paths', counts(0, 'paths', extension='completed'),
             'fields|a plain count'),
            ('supported_paths', at_least(0), 'for a search completed'),
            ('supported_paths', counts(0, 'supported_paths', level='label_intersection'),
             'the index\'s best support'),
            ('supported_paths_label', counts(0, 'supported_paths', level='record_verified'),
             'the level asked for'),
            ('supported_paths', counts(0, 'paths', value=1), 'every supported path is a complete'),
            ('supported_paths_predicate_pruned', counts(0, 'paths', relation='exact'),
             'walks a predicate pruned'),
            ('supported_paths_predicate_either',
             both(at_least(0), counts(0, 'supported_paths', branches_pruned_by_predicate=1)),
             'does not prune'),
            ('supported_paths', work(0, anchor_rows=100), 'rows read'),
            ('supported_paths_predicate_either', work(0, mirror_rows=3), 'one strand'),
            ('supported_paths_predicate_mirror', work(0, predicate_rows=11),
             'reads no row of its own'),
            ('supported_paths_work_budget', counts(0, 'supported_paths', search='completed'),
             'for a search completed|the supported paths at_least'),
            ('supported_paths_partial', entry(0, cut={'reason': 'max_annotation_work'}),
             'a stop at max_annotation_work'),
            ('paths_max_paths_partial', entry(0, cut={'reason': 'max_annotation_work'}),
             'max_annotation_work cuts the supported paths'),
            ('supported_paths_work_budget', entry(0, withheld={'reason': 'count_above_threshold'}),
             'the exact supported paths above max_paths'),
            ('supported_paths_work_budget', entry(0, withheld={'reason': 'output_budget'}),
             'did not fit the account'),
            ('supported_paths_memory',
             lambda a: a['patterns'][5]['rows_refused'][0].update(phase='discovery'), 'not in'),
            ('supported_paths_predicate_pruned',
             results(0, lambda r: r.update(selection_strands=['reverse_complement'])),
             'the own row only'),
            ('supported_paths_predicate_mirror',
             results(0, lambda r: r.update(selection_strands=['context'])), 'on the own row'),
            ('supported_paths', entry(0, notes=['paths_later_increment']),
             'paths_later_increment exactly'),
            ('supported_paths_count', entry(0, notes=[]), 'projection_not_read exactly'),
            ('supported_paths_work_budget', entry(0, notes=['projection_not_read']),
             'projection_not_read exactly'),
            ('supported_paths_label', entry(0, notes=[]), 'label_intersection_only exactly'),
            ('supported_paths_work_budget', entry(0, notes=['low_complexity_pattern']),
             'beside a budget stop'),
            ('supported_paths', lambda a: a['patterns'][0]['timing'].update(
                label_discovery_ms=0.1), 'fields'),
            ('supported_paths', lambda a: a['limits'].pop('supported_paths_level'), 'fields'),
            ('supported_paths', lambda a: a['limits'].update(supported_paths_level='best2'),
             'as requested'),
            ('supported_paths_predicate_label',
             lambda a: a['patterns'][0]['selection'].update(support='record_verified'),
             'searched at'),
            ('supported_paths_predicate_pruned', counts(0, 'tested', value=2),
             'every raw context tested|selected <= tested'),
            ('supported_paths_label', entry(0, placement='record'), 'fields|counts.occurrences'),
            ('supported_paths_count', lambda a: a['patterns'][0]['work'].pop('anchor_rows'),
             'fields'),
            ('supported_paths_count', lambda a: a['patterns'][0].pop('rows_refused'), 'fields'),
            ('supported_paths_memory', entry(4, stop=None), 'a max_memory cut states its stop'),
            ('supported_paths_above_threshold',
             entry(0, withheld={'reason': 'selected_above_threshold'}),
             'a predicate\'s withheld or cut without a predicate'),
            ('supported_paths_partial', entry(0, cut={'reason': 'max_predicate_work'}),
             'a predicate\'s withheld or cut without a predicate'),
        ]
        self.assertRefusesMutations(cases)
        # requests the route refuses (SPEC §20.2, §25.2) are never answered 200
        request, _ = self.bodies['supported_paths']
        for change, says in (({'supported_paths_level': 'label_intersection',
                               'require_support': 'record_verified'}, 'record_verified'),
                             ({'supported_paths_level': 'all'}, 'supported_paths_level'),
                             ({'predicate_scope': 'motif'}, 'needs a predicate')):
            with self.subTest(change=change):
                with self.assertRaisesRegex(AssertionError, says):
                    self.request_200(Checker(Raising(), 'request'), {**request, **change})

    def test_motif_rules_refuse_what_v1_never_answers(self):
        """SPEC §25: the stored bodies pass, and each rule refuses an answer that breaks it -- a
        value that is not the normal form on the union, an absence claim from an incomplete
        pass, a tested_contexts value the labels found do not decide (or an undecided one they
        do), every_context where a context went untested, for a pattern longer than k with
        anchors or for one without an instance, no_instance for a pattern with contexts (and a
        pattern without an instance decided any other way, or undecided), the labels found out
        of label order, with strands their rows do not give or more contexts than tested, the
        motif's units outside the selection's, the block or the scope missing or where no motif
        was asked, and the capabilities' scopes."""
        def motif(i, **fields):
            return lambda a: a['patterns'][i]['motif'].update(fields)

        def label(i, j, **fields):
            return lambda a: a['patterns'][i]['motif']['labels_present'][j].update(fields)

        names = [n for n in self.fixtures if n.startswith('motif_')]
        self.assertGreaterEqual(len(names), 4)
        for name in names:
            self.check_stored(name)
        cases = [
            ('motif_context', motif(0, selected=False), 'the normal form on the union'),
            ('motif_context', motif(0, labels_absent=0), 'not found'),
            ('motif_tested_contexts', motif(0, labels_absent=1), 'no absence claim'),
            ('motif_tested_contexts', motif(0, decided_by='every_context'),
             'exact with every context|every context tested'),
            ('motif_tested_contexts', motif(0, selected=False), 'definite on the labels found'),
            ('motif_tested_contexts', motif(0, selected=None, decided_by=None),
             'undecided although'),
            ('motif_either', label(0, 0, strands='either'), 'one row serves both'),
            ('motif_context', label(0, 0, strands='both'), 'the own row only'),
            ('motif_either', lambda a: a['patterns'][0]['motif']['labels_present'].reverse(),
             'label order'),
            ('motif_context', label(0, 0, contexts={'value': 13, 'relation': 'exact',
                                                    'unit': 'graph_contexts'}),
             'at most the contexts tested'),
            ('motif_long_patterns', motif(0, selected=True, decided_by='every_context',
                                          untested=None, labels_absent=1),
             'the pass completed on at least one context'),
            # the pattern without an instance: its value the normal form on the empty set, every
            # label absent, nothing untested; never every_context, constant or undecided
            ('motif_long_patterns', motif(1, selected=False), 'the normal form on the empty set'),
            ('motif_long_patterns', motif(1, labels_absent=0), 'every label of the normal form'),
            ('motif_long_patterns', motif(1, untested='not_started'),
             'nothing to test, nothing found'),
            ('motif_long_patterns', motif(1, decided_by='every_context'),
             'the pass completed on at least one context'),
            ('motif_long_patterns', motif(1, decided_by='constant'), 'a constant normal form'),
            ('motif_long_patterns', motif(1, selected=None, decided_by=None,
                                          untested='not_started', labels_absent=None),
             'a pattern without an instance says so'),
            # no_instance only without an instance
            ('motif_context', motif(0, decided_by='no_instance'), 'the raw count exact 0'),
            ('motif_context', lambda a: a['patterns'][0]['work'].update(motif_units=10 ** 9),
             'part of predicate_units'),
            ('motif_context', lambda a: a['patterns'][0].pop('motif'), 'fields'),
            ('motif_context', lambda a: a['predicate'].update(motif_scope='shard_context'),
             'not in'),
            ('motif_context', lambda a: a['predicate'].pop('motif_scope'), 'fields'),
            ('predicate_context', lambda a: a['predicate'].update(motif_scope='shard_motif'),
             'fields'),
            ('motif_context', lambda a: a['limits'].update(predicate_scope='context'),
             'as requested'),
        ]
        self.assertRefusesMutations(cases)
        b = copy.deepcopy(self.capabilities_of('masked'))
        b['predicate']['scopes'] = ['context']
        with self.assertRaisesRegex(AssertionError, 'scopes'):
            Checker(Raising(), 'capabilities').block(b, 'pattern', False)

    def test_documents_state_the_built_deadline(self):
        """The route's default deadline and its cap, as the SPEC and the service's request
        (PROMPT-search-service-pattern.md §1 and §3.2) state them, are the capabilities'."""
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
                                ('--pattern-max-time-ms', caps['caps']['time_budget_ms']),
                                # the two other defaults a host sets apart from the caps
                                ('--pattern-default-max-anchors', caps['default_max_anchors']),
                                ('--pattern-max-anchors', caps['caps']['max_anchors']),
                                ('--pattern-default-max-labels', caps['default_max_labels']),
                                ('--pattern-max-labels', caps['caps']['max_labels'])):
                self.assertIn(f'| `{flag}` | {value:,} |', spec, flag)
        # a default is at most its cap (the server refuses a larger one at start-up)
        self.assertLessEqual(caps['default_max_anchors'], caps['caps']['max_anchors'])
        self.assertLessEqual(caps['default_max_labels'], caps['caps']['max_labels'])

    def test_capabilities_documents_keep_a_kibibyte(self):
        """Every capabilities document of every fixture server -- GET /traverse/capabilities is
        the one a service's MCP tool returns in one piece, under CAPABILITIES_MAX_BYTES -- keeps
        CAPABILITIES_ROOM bytes under that ceiling as the server writes it (served_bytes, from
        above), so that later additions fit. The paths a multi-graph server names are measured
        as stored ({work}/...): a host's own paths are its own."""
        if os.path.isfile(MCP_TOOLS):
            with open(MCP_TOOLS, encoding='utf-8') as f:
                m = re.search(r'^CAPABILITIES_MAX_BYTES = (\d+) \* 1024$', f.read(), re.M)
            self.assertTrue(m, 'CAPABILITIES_MAX_BYTES in mcp_tools.py')
            self.assertEqual(CAPABILITIES_MAX_BYTES, int(m.group(1)) * 1024)
        budget = CAPABILITIES_MAX_BYTES - CAPABILITIES_ROOM
        servers = set()
        for name, f in self.fixtures.items():
            if f['method'] != 'GET':
                continue
            servers.add(f['server'])
            size = served_bytes(self.bodies[name][1])
            with self.subTest(fixture=name):
                self.assertLessEqual(size, budget, f'{name}: {size} bytes, {budget - size} '
                                                   f'over the budget of {budget}')
        # every server a fixture runs on but the ones answering no GET (variants of a measured
        # server by a flag: masked_small_predicate_cap's document is the masked one's with a
        # smaller caps.max_predicate_labels, one digit shorter; multi_mixed's is the multi one's
        # with the hash pair's graph_summary entry, the hash server's reason in it)
        self.assertEqual({f['server'] for f in self.fixtures.values()}
                         - {'masked_no_map', 'unmasked_unchecked', 'masked_small_predicate_cap',
                            'multi_mixed'},
                         servers)
        # the measure: exact on what it does not bound from above
        self.assertEqual(len(json.dumps({'a': [1, 'b', None, True, {}]}, separators=(',', ':'))),
                         served_bytes({'a': [1, 'b', None, True, {}]}))
        self.assertEqual(len('{"x":"\\u00a7"}'), served_bytes({'x': '§'}))
        self.assertEqual(len('{"f":}') + 24, served_bytes({'f': 0.5}))


if __name__ == '__main__':
    unittest.main()
