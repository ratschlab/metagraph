#!/usr/bin/env python3
"""Build (or --check) the POST /pattern fixture bodies of docs/SPEC-pattern-search.md:
api/python/tests/data/traverse/pattern/<name>/{request.json,answer.json}, with index.json
(what each fixture is: server, method, path, HTTP status) and README.md (one line each).

The bodies are what a server_query of the given binary answers on copies of the mini index
(build/mini_refseq, scripts/traversal/build_mini_refseq.sh), so that the search service can
test its pattern job against the frozen contract (pattern_contract_version 1) without a host
(PROMPT-search-service-pattern.md §2, "fixtures"). The mini index is only read, never
modified: everything is built in a work directory.

Servers (each a server_query on 127.0.0.1, started and stopped by this script):

  masked     the mini index with its dummy-edge mask: a copy of the mini's graph given its
             .edgemask by `metagraph transform --mask-dummy` (DESIGN-pattern-search.md §4, the
             one-time step for a host; the .dbg is unchanged, checked, so the mini's annotation
             and .seqs serve it as they are), as integration_tests/test_pattern.py makes it
  unmasked   the mini index as built (no .edgemask): the route answers mask_required
  built_at_load
             the mini index as built, served with --pattern-build-mask: the same mask built in
             memory at start-up (capabilities mask: built_at_load)
  multi      a multi-graph server (graphs CSV) carrying the masked copy as `mini_refseq`
  primary    a PRIMARY index of the mini's NDM-carrying records (562.fa, 573.fa), served
             wrapped in CanonicalDBG: orientations instead of strands, no suffix scope; its
             column annotation has no budget-aware decode (annotation: unbudgeted)
  masked_no_map
             the masked copy served with --no-coord-mapping: coordinates without the record
             mapping (output.labels "all" places nothing: placement global)
  hash       a hash graph (build --graph hash) of one mini record (1296536.fa) with a column
             annotation: a graph the engine does not recognise (representation_unsupported;
             review of 2026-10-07, C2-01)

Two fixtures are HAND-MADE (index.json "hand_made": true): the 503 bodies a server cannot be
made to produce on demand (an answer that overran its finalisation reserve; a request during
the index load). Their text is the code's (src/cli/pattern.cpp PatternDelivery::check,
server_utils.cpp process_request) and the unit test checks it against the source. Two refusal
codes have no fixture at all: primary_unwrapped (server_query and the CLI always wrap a PRIMARY
graph in CanonicalDBG) and alphabet_unsupported (a graph of another alphabet does not load in
this build).

Volatile values: every answer is stored as the server wrote it, but `timing.elapsed_ms` (and
every other `timing` value of a pattern), the per-process `server_instance`, and the counts and
work of the FIRST `determinism: time_limited` pattern of an answer (the one the clock stopped:
how far it got; when the clock cut its release, also its `returned` and `results`) depend on the
run; --check blanks them on both sides, and a regeneration keeps
a file whose blanked content did not change (so that rerunning does not churn the files). The
patterns after it are stopped by the same budget (SPEC §7.6: unknown counts, work zero), the
same in every run, and compared as they are (review of 2026-10-07, C2-02). Paths
under the work directory are written as {work}/... (they would otherwise name a temporary
directory). Each fixture also asserts the situation it exists to show (`expect` below), so a
server change cannot silently turn a fixture into a different case.

usage:
    pattern_fixtures.py [--metagraph BIN] [--mini DIR] [--dir DIR] [--work DIR] [--check]
                        [--only NAME ...]
  --metagraph  the binary (default: build/metagraph of this checkout)
  --mini       the mini index (default: build/mini_refseq)
  --work       keep the built copies there (reused when present) instead of a temp dir
  --check      write nothing; exit 1 when a regenerated file differs (request.json,
               README.md and index.json byte for byte, answer.json after blanking)
"""

import argparse
import copy
import filecmp
import glob
import json
import os
import shutil
import socket
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.request

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DEFAULT_DIR = os.path.join(REPO, 'api', 'python', 'tests', 'data', 'traverse', 'pattern')
DEFAULT_MINI = os.path.join(REPO, 'build', 'mini_refseq')
DEFAULT_BIN = os.path.join(REPO, 'build', 'metagraph')

MINI_K = 31
MINI_GRAPH = 'graph_k31.dbg'
MINI_ANNO = 'annotation.relaxed.relabeled.row_diff_brwt_coord.annodbg'
# served beside the graph and the annotation: the record mapping and the row-diff sidecars
MINI_SIDECARS = ['annotation.relaxed.relabeled.seqs', 'graph_k31.dbg.anchors',
                 'graph_k31.dbg.rd_succ']
PRIMARY_RECORDS = ['562.fa', '573.fa']

INDEX_NAME = 'mini_refseq'
INDEX_RELEASE = 'fixtures-mini_refseq'
PRIMARY_NAME = 'mini_primary'
MULTI_GRAPH_NAME = 'mini_refseq'

# ----------------------------------------------------------------------- the patterns
#
# Real primers where one exists, so that a reader sees what the route is for:
# blaNDM-1 (carried by the mini's ndm_* records on both strands), the 16S V4 primers (in every
# genome of the mini) and a palindromic EcoRI-centred 12-mer.
NDM_F = 'GGTTTGGCGATCTGGTTTTC'       # blaNDM-1 forward primer, 20 nt, 40 bits
NDM_R = 'CGGAATGGCTCATCACGATC'       # blaNDM-1 reverse primer
NDM_F11 = NDM_F[:11]                 # 22 bits: below the default floor of 24 outside suffix
NDM_40 = 'ATGGAATTGCCCAATATTATGCACCCGGTCGCGAAGCTGA'   # the first 40 nt of blaNDM-1 (L > k)
P515F = 'GTGYCAGCMGCCGCGGTAA'        # 16S V4 515F (Parada), IUPAC, 36 bits
P806R = 'GGACTACNVGGGTWTCTAAT'       # 16S V4 806R (Apprill), IUPAC with an N
PAL12 = 'GCCGAATTCGGC'               # palindromic (P = rc(P)), the EcoRI site at its centre
LOW = 'GCCGCCGCCGCCGCCGCCGC'         # a repeat sdust flags (low_complexity_pattern)
# random 20- and 40-mers absent from the mini (asserted: exact zeros)
ABSENT_20 = 'CTAGGAGATGGGCCAGCTAC'
ABSENT_40 = 'GATAGAGAACTCGAGAGAGGTTCCACCTTCATATTGAATT'
# 12 bases of NDM-F's first 13 around a run of 19 Ns: 24 bits, admitted, but ~2.5e7 range
# steps on the mini (about 1.5 s): stopped by any budget of a few milliseconds. (It was 12 Ns
# then NDM-F's first 12 bases, until the engine learnt to skip a pattern's leading and
# trailing N runs (review of 2026-10-07, X-EFFICIENCY-01): that one now takes 125 steps)
HEAVY = NDM_F[:2] + 'N' * 19 + NDM_F[3:13]
# a 10-mer whose suffix-scope count needs a mask scan (the prefix of a record's first k-mer):
# discovery takes 20 range steps, its scan one more, so max_steps 20 stops in the scan with
# bounds [lower, upper] (review of 2026-10-07, X-TESTS-03: no body carried bounds or mask_scan)
SCAN_10 = 'ATGCCGGTGA'
SCAN_STEPS = 20
# the record the hash graph is built from (a graph the engine does not recognise)
HASH_RECORD = '1296536.fa'

# the server's defaults (--pattern-* flags, DESIGN-pattern-search.md §5.3), for the hand-made
# 503 and the expectations
FINALIZE_MS = 250
TINY_BUDGET_MS = FINALIZE_MS + 1
# a step stop whose release then meets the work time (review of 2026-10-07, C1-03): ACG in
# suffix scope (an exact pattern of L <= k there: exempt from the floor) has ~1.2e5 contexts
# on the mini, 3 range steps stop its discovery, and 150 ms of work time are always spent in
# its release before max_contexts (10,000) are out: the estimated time to write the results
# built (SPEC §7.6, ~0.15 ms per KB at the default rates) alone passes it after ~5,500 of them
STEP_THEN_TIME_STEPS = 3
STEP_THEN_TIME_BUDGET_MS = FINALIZE_MS + 150


# ----------------------------------------------------------------------- expectations

def counts_of(entry):
    c = entry['counts']
    return c.get('contexts') or c['anchors']


def check(condition, what):
    if not condition:
        raise AssertionError(what)


def exact(value):
    def run(entry):
        c = counts_of(entry)
        check((c['relation'], c['value']) == ('exact', value), (c['relation'], c['value']))
    return run


def withheld(reason):
    def run(entry):
        check(entry['withheld'] == {'reason': reason}, entry['withheld'])
    return run


def cut(reason):
    def run(entry):
        check(entry['cut'] == {'reason': reason}, entry['cut'])
    return run


def field(name, value):
    def run(entry):
        check(entry.get(name) == value, (name, entry.get(name)))
    return run


def counted(name, relation, value=None):
    """counts.<name> (labels, occurrences) with |relation| (and |value| when given)."""
    def run(entry):
        c = entry['counts'][name]
        check(c['relation'] == relation and (value is None or c['value'] == value), (name, c))
    return run


def labelled(entry):
    """output.labels "all": every result carries its labels read in full."""
    for r in entry['results']:
        check(r['labels_status'] == 'complete' and isinstance(r['labels'], list), r)


def slot_error(code):
    def run(entry):
        check(entry['error']['code'] == code, entry.get('error'))
    return run


def complete(entry):
    check(entry['retrieval_complete'] is True and entry['withheld'] is None
          and entry['cut'] is None and entry['stop'] is None, entry)


def unknown_counts(entry):
    """Every graph count of the entry unknown and null: the total, suffix, every offset and
    strand (SPEC §7.6: a pattern after a budget stop)."""
    c = counts_of(entry)
    parts = [c, c.get('suffix')] + list(c.get('by_offset', {}).values()) \
        + list(c.get('by_strand', c.get('by_orientation', {})).values())
    for x in parts:
        if x is not None:
            check(x['relation'] == 'unknown' and x['value'] is None, x)


def work_zero(entry):
    check(all(v == 0 for v in entry['work'].values()), entry['work'])


def relation(value):
    """counts.contexts (or anchors) with this relation."""
    def run(entry):
        check(counts_of(entry)['relation'] == value, counts_of(entry))
    return run


def not_exact(entry):
    """A stop in discovery leaves no exact total (SPEC §7.4); nothing of the annotation is read."""
    check(counts_of(entry)['relation'] != 'exact', counts_of(entry))
    for name in ('labels', 'occurrences'):
        check(entry['counts'][name]['relation'] == 'unknown', entry['counts'][name])


def bounded(entry):
    """relation bounds with value == lower <= upper."""
    c = counts_of(entry)
    check(c['relation'] == 'bounds' and c['value'] == c['lower'] <= c['upper'], c)


def refused(code):
    def run(answer):
        check(set(answer) == {'error', 'code'} and answer['code'] == code, answer)
    return run


def expect_all(*checks):
    def run(value):
        for c in checks:
            c(value)
    return run


def entries(*per_entry):
    """One check per pattern entry, in request order."""
    def run(answer):
        check(len(answer['patterns']) == len(per_entry), len(answer['patterns']))
        for entry, c in zip(answer['patterns'], per_entry):
            c(entry)
    return run


def caps_block(available, reason=None, mask=None):
    def run(doc):
        b = doc['pattern']
        check(b['pattern_contract_version'] == 1, b)
        check(b['available'] is available and b['unavailable_reason'] == reason, b)
        if mask is not None:
            check(b['mask'] == mask, b)
    return run


def only_k(doc):
    """A graph the engine does not recognise: graph_mode null, only k set (SPEC §10.2)."""
    b = doc['pattern']
    check(b['graph_mode'] is None and isinstance(b['k'], int), b)
    for f in ('alphabet', 'strand_stated', 'mask', 'scopes', 'placement', 'support', 'annotation'):
        check(b[f] is None, (f, b[f]))


def features(listed):
    def run(doc):
        check(('pattern' in doc['features']) is listed, doc['features'])
        check(('pattern' in doc['routes']) is listed, doc['routes'])
    return run


# ----------------------------------------------------------------------- the fixtures

class Fixture:
    """One request to one server, the HTTP status it must get, and what its answer shows."""

    def __init__(self, name, server, method, path, body, status, shows, expect=None,
                 hand_made=None, headers=None):
        self.name = name
        self.server = server
        self.method = method
        self.path = path
        self.body = body
        self.status = status
        self.shows = shows
        self.expect = expect
        # the stored answer of a hand-made fixture (no server produces it on demand)
        self.hand_made = hand_made
        # response headers a client acts on (Retry-After of the 503 during the load)
        self.headers = headers or {}


def post(name, server, body, status, shows, expect=None, **kw):
    return Fixture(name, server, 'POST', '/pattern', body, status, shows, expect, **kw)


def get(name, server, path, shows, expect=None):
    return Fixture(name, server, 'GET', path, None, 200, shows, expect)


def p(text, kind='dna', ident=None):
    """One request pattern."""
    out = {kind: text}
    if ident is not None:
        out['id'] = ident
    return out


FIXTURES = [
    # ---------------------------------------------------------------- capabilities
    get('capabilities', 'masked', '/capabilities',
        'GET /capabilities of a single-graph server whose graph has its mask: `pattern` in '
        'features and routes, the block available (basic, mask file, placement record)',
        expect_all(features(True), caps_block(True, mask='file'))),
    get('traverse_capabilities', 'masked', '/traverse/capabilities',
        'GET /traverse/capabilities (the document the service probe reads) on the same server: '
        'the same `pattern` block',
        caps_block(True, mask='file')),
    get('capabilities_mask_absent', 'unmasked', '/capabilities',
        'the mini index as built, without its .edgemask: the feature and route listed, the '
        'block available false, unavailable_reason mask_required, mask absent',
        expect_all(features(True), caps_block(False, 'mask_required', 'absent'))),
    get('traverse_capabilities_mask_absent', 'unmasked', '/traverse/capabilities',
        'the same unmasked server on the probe route',
        caps_block(False, 'mask_required', 'absent')),
    get('capabilities_built_at_load', 'built_at_load', '/capabilities',
        'the mini index as built, served with --pattern-build-mask: the block available, mask '
        'built_at_load (the mask built in memory at start-up), otherwise as with the file',
        expect_all(features(True), caps_block(True, mask='built_at_load'))),
    get('traverse_capabilities_built_at_load', 'built_at_load', '/traverse/capabilities',
        'the same server on the probe route',
        caps_block(True, mask='built_at_load')),
    get('capabilities_multi_graph', 'multi', '/capabilities',
        'a multi-graph server: no `pattern` feature or route; the block says '
        'multi_graph_later_increment and nothing else',
        expect_all(features(False), caps_block(False, 'multi_graph_later_increment'))),
    get('traverse_capabilities_multi_graph', 'multi',
        '/traverse/capabilities?graph=' + MULTI_GRAPH_NAME,
        'the multi-graph server probed for one graph: the same reduced block',
        caps_block(False, 'multi_graph_later_increment')),
    get('traverse_capabilities_primary', 'primary', '/traverse/capabilities',
        'a PRIMARY index (wrapped in CanonicalDBG): graph_mode primary, scopes [any_offset], '
        'strand_stated false, placement none_canonical, annotation unbudgeted (column)',
        caps_block(True, mask='file')),
    get('capabilities_representation_unsupported', 'hash', '/capabilities',
        'a graph the engine does not recognise (a hash graph): the feature and route listed, '
        'the block available false, unavailable_reason representation_unsupported, graph_mode '
        'null and only k set',
        expect_all(features(True), caps_block(False, 'representation_unsupported'), only_k)),
    get('traverse_capabilities_representation_unsupported', 'hash', '/traverse/capabilities',
        'the same hash-graph server on the probe route',
        expect_all(caps_block(False, 'representation_unsupported'), only_k)),

    # ---------------------------------------------------------------- count
    post('count', 'masked',
         {'patterns': [p(NDM_F, ident='NDM-F'), p(NDM_R, ident='NDM-R'), p(ABSENT_20, ident='absent')],
          'mode': 'count'}, 200,
         'mode count, scope any_offset and strands both by default: exact context counts with '
         'suffix, by_offset (every offset, zeros included) and by_strand; labels and '
         'occurrences unknown; no results; an absent primer counts exact 0',
         entries(exact(24), exact(24), exact(0))),
    post('count_suffix', 'masked',
         {'patterns': [p(NDM_F), p(NDM_F11)], 'mode': 'count', 'scope': 'suffix'}, 200,
         'scope suffix: only offset k - L, absence_scope suffix_only; an exact 11-mer (22 '
         'bits) is admitted in this scope however short',
         entries(exact(2), exact(8))),
    post('count_forward', 'masked',
         {'patterns': [p(NDM_F)], 'mode': 'count', 'strands': 'forward'}, 200,
         'strands forward: P only, strands ["+"], by_strand {"+"}',
         entries(expect_all(exact(12), field('strands', ['+'])))),
    post('count_reverse', 'masked',
         {'patterns': [p(NDM_F)], 'mode': 'count', 'strands': 'reverse'}, 200,
         'strands reverse: rc(P) only, strands ["-"]',
         entries(expect_all(exact(12), field('strands', ['-'])))),
    post('count_iupac', 'masked',
         {'patterns': [p(P515F, 'iupac', '515F'), p(P806R, 'iupac', '806R')], 'mode': 'count'},
         200,
         'IUPAC patterns (the 16S V4 primers 515F and 806R): kind iupac, fractional '
         'information_bits',
         entries(exact(26), exact(24))),
    post('count_palindrome', 'masked',
         {'patterns': [p(PAL12)], 'mode': 'count'}, 200,
         'a palindromic pattern: searched once, palindromic true, strands ["="], by_strand '
         '{"both"}',
         entries(expect_all(exact(220), field('palindromic', True), field('strands', ['='])))),
    post('count_low_complexity', 'masked',
         {'patterns': [p(LOW)], 'mode': 'count'}, 200,
         'a repeat sdust flags: notes ["low_complexity_pattern"] (here with no context)',
         entries(field('notes', ['low_complexity_pattern']))),
    post('count_long', 'masked',
         {'patterns': [p(NDM_40), p(ABSENT_40)], 'mode': 'count'}, 200,
         'L > k in mode count: scope long, counts.anchors (anchor window [0, k)) and '
         'counts.paths unknown, anchor_information_bits, note paths_later_increment; no '
         'anchor: paths exact 0',
         entries(expect_all(exact(2), field('scope', 'long')), exact(0))),
    post('information_floor', 'masked',
         {'patterns': [p(NDM_F11, ident='short'), p('NNNNNNNNACGTACGTAC', 'iupac', 'degenerate'),
                       p('N' * 20 + NDM_F, 'iupac', 'long_window'), p(NDM_F, ident='ok')],
          'mode': 'count'}, 200,
         'error slots information_below_floor (22 bits; 20 bits; a long pattern whose anchor '
         'window has 22): the slot keeps the pattern\'s description, the others are answered',
         entries(*[slot_error('information_below_floor')] * 3, exact(24))),
    post('bad_alphabet', 'masked',
         {'patterns': [p('ACGUACGTACGTACGTACGT', ident='u'), p('ACGRACGTACGTACGTACGT', ident='r'),
                       p('ACGTXACGTACGTACGTACG', 'iupac', 'x'), p('', ident='empty'),
                       p(NDM_F, ident='ok')],
          'mode': 'count'}, 200,
         'error slots bad_alphabet (U, an IUPAC code in a dna pattern, X, an empty pattern): '
         'only id, kind and error; the last pattern is answered',
         entries(*[slot_error('bad_alphabet')] * 4, exact(24))),
    post('clamped', 'masked',
         {'patterns': [p(NDM_F)], 'mode': 'count', 'max_contexts': 1000000,
          'max_steps': 1000000000, 'time_budget_ms': 1000000}, 200,
         'request values above the server caps are lowered and listed in limits.clamped '
         '{field, requested, effective}, never refused',
         lambda a: check(len(a['limits']['clamped']) == 3, a['limits'])),

    # ---------------------------------------------------------------- retrieval
    post('all_or_count', 'masked',
         {'patterns': [p(NDM_F, ident='NDM-F'), p(ABSENT_20, ident='absent')]}, 200,
         'mode and output omitted: all_or_count with output.labels "none"; every context '
         'returned (kmer, instance, offset, strand, node, row) in (node, offset, orientation) '
         'order, retrieval_complete true; the absent primer: results [] and complete',
         entries(expect_all(exact(24), complete), expect_all(exact(0), complete))),
    post('all_or_count_iupac', 'masked',
         {'patterns': [p(P515F, 'iupac', '515F')], 'mode': 'all_or_count',
          'output': {'labels': 'none'}}, 200,
         'an IUPAC pattern retrieved: each result\'s instance is the matched variant (Y, M '
         'resolved), rc(P)\'s instances on strand "-"',
         entries(expect_all(exact(26), complete))),
    post('all_or_count_suffix', 'masked',
         {'patterns': [p(NDM_F), p(PAL12)], 'scope': 'suffix',
          'output': {'labels': 'none'}}, 200,
         'retrieval in scope suffix (offset k - L only, absence_scope suffix_only); the '
         'palindrome\'s contexts carry strand "="',
         entries(expect_all(exact(2), complete), expect_all(exact(11), complete))),
    post('all_or_count_withheld', 'masked',
         {'patterns': [p(NDM_F)], 'max_contexts': 10, 'output': {'labels': 'none'}}, 200,
         'the exact count (24) above max_contexts (10): results withheld, '
         'withheld.reason count_above_threshold, returned 0, no cut',
         entries(expect_all(exact(24), withheld('count_above_threshold'), field('returned', 0)))),
    post('partial_cut', 'masked',
         {'patterns': [p(NDM_F)], 'mode': 'partial', 'max_contexts': 5,
          'output': {'labels': 'none'}}, 200,
         'mode partial: the first 5 contexts in answer order, cut.reason max_contexts, '
         'retrieval_complete false, the count exact',
         entries(expect_all(exact(24), cut('max_contexts'), field('returned', 5)))),
    post('stop_at_threshold', 'masked',
         {'patterns': [p(NDM_F)], 'max_contexts': 5, 'stop_at_threshold': True,
          'output': {'labels': 'none'}}, 200,
         'stop_at_threshold: discovery stops once the count passes max_contexts: at_least, '
         'stop {discovery, max_contexts}, withheld threshold_crossed',
         entries(expect_all(withheld('threshold_crossed'),
                            field('stop', {'phase': 'discovery', 'reason': 'max_contexts'})))),
    post('max_steps', 'masked',
         {'patterns': [p(NDM_F), p(NDM_R)], 'max_steps': 50, 'output': {'labels': 'none'}},
         200,
         'max_steps reached in discovery: at_least counts, stop {discovery, max_steps}, '
         'withheld discovery_budget; the request\'s budget is spent, so the next pattern '
         'answers unknown counts with the same stop',
         entries(withheld('discovery_budget'), withheld('discovery_budget'))),
    post('max_steps_partial', 'masked',
         {'patterns': [p(NDM_F), p(NDM_R)], 'mode': 'partial', 'max_steps': 50,
          'output': {'labels': 'none'}}, 200,
         'the same stop in mode partial: what discovery found is delivered, cut.reason '
         'max_steps; the second pattern returns nothing, cut max_steps',
         entries(cut('max_steps'), expect_all(cut('max_steps'), field('returned', 0)))),
    post('max_steps_bounds', 'masked',
         {'patterns': [p(SCAN_10, ident='scan')], 'mode': 'count', 'scope': 'suffix',
          'max_steps': SCAN_STEPS}, 200,
         'max_steps reached in a mask scan (scope suffix, a record-start pattern whose count '
         'needs one): stop {mask_scan, max_steps}, the count bounds {value = lower, upper} with '
         'the strand that finished exact; bounds + exact = bounds',
         entries(expect_all(bounded, field('stop', {'phase': 'mask_scan',
                                                   'reason': 'max_steps'})))),
    post('max_steps_bounds_withheld', 'masked',
         {'patterns': [p(SCAN_10, ident='scan')], 'mode': 'all_or_count', 'scope': 'suffix',
          'max_steps': SCAN_STEPS, 'output': {'labels': 'none'}}, 200,
         'the same stop in all_or_count: bounds, withheld discovery_budget',
         entries(expect_all(bounded, withheld('discovery_budget'),
                            field('stop', {'phase': 'mask_scan', 'reason': 'max_steps'})))),
    post('deadline', 'masked',
         {'patterns': [p(HEAVY, 'iupac', 'heavy'), p(NDM_F, ident='NDM-F')],
          'time_budget_ms': TINY_BUDGET_MS, 'output': {'labels': 'none'}}, 200,
         'a 251 ms budget (1 ms of work before the 250 ms finalisation reserve): stop '
         '{discovery, time}, withheld deadline, determinism time_limited; the next pattern '
         'unknown (counts and work of time_limited patterns vary between runs)',
         entries(expect_all(withheld('deadline'), field('determinism', 'time_limited'),
                             field('stop', {'phase': 'discovery', 'reason': 'time'}),
                             not_exact),
                 # the budget is the request's: the next pattern unknown, work zero, the same
                 # stop (SPEC §7.6; compared as stored, not blanked)
                 expect_all(withheld('deadline'), field('determinism', 'time_limited'),
                            field('stop', {'phase': 'discovery', 'reason': 'time'}),
                            unknown_counts, work_zero))),
    post('deadline_partial', 'masked',
         {'patterns': [p(HEAVY, 'iupac', 'heavy')], 'mode': 'partial',
          'time_budget_ms': TINY_BUDGET_MS, 'output': {'labels': 'none'}}, 200,
         'the deadline in mode partial: nothing is released after a time stop in discovery '
         '(its membership would depend on the machine): returned 0, cut.reason time',
         entries(expect_all(cut('time'), field('returned', 0), not_exact,
                            field('stop', {'phase': 'discovery', 'reason': 'time'})))),
    post('max_steps_then_time', 'masked',
         {'patterns': [p('ACG', ident='steps'), p('GGA', ident='after')], 'mode': 'partial',
          'scope': 'suffix', 'max_steps': STEP_THEN_TIME_STEPS,
          'time_budget_ms': STEP_THEN_TIME_BUDGET_MS, 'output': {'labels': 'none'}}, 200,
         'partial, discovery stopped by max_steps, then the release meets the work time: stop '
         'keeps the first stop {discovery, max_steps}, the time shows only as cut.reason time '
         'and determinism time_limited (how many were released varies between runs); the next '
         'pattern answers the sticky max_steps stop, its empty release past the work time: cut '
         'time, returned 0, time_limited',
         entries(expect_all(field('stop', {'phase': 'discovery', 'reason': 'max_steps'}),
                            cut('time'), field('determinism', 'time_limited'), not_exact),
                 expect_all(field('stop', {'phase': 'discovery', 'reason': 'max_steps'}),
                            cut('time'), field('returned', 0),
                            field('determinism', 'time_limited'), unknown_counts,
                            work_zero))),
    post('long', 'masked',
         {'patterns': [p(NDM_40, ident='NDM-40'), p(ABSENT_40, ident='absent')],
          'output': {'labels': 'none'}}, 200,
         'L > k in all_or_count: anchors counted, paths unknown, withheld '
         'paths_later_increment; with no anchor the empty answer is complete (anchors and '
         'paths exact 0)',
         entries(expect_all(exact(2), withheld('paths_later_increment')),
                 expect_all(exact(0), complete))),
    post('primary_any_offset', 'primary',
         {'patterns': [p(NDM_F, ident='NDM-F')], 'output': {'labels': 'none'}}, 200,
         'a PRIMARY index: results carry orientation (forward, reverse) instead of strand, '
         'by_orientation instead of by_strand, note strand_unknown_canonical; node is the '
         'wrapper id',
         entries(expect_all(exact(24), complete, field('strands', ['forward', 'reverse'])))),
    post('primary_suffix', 'primary',
         {'patterns': [p(NDM_F)], 'mode': 'count', 'scope': 'suffix'}, 200,
         'scope suffix on a PRIMARY index: error slot scope_unsupported (a virtual suffix is a '
         'stored prefix); any_offset is complete there',
         entries(slot_error('scope_unsupported'))),

    # ---------------------------------------------------------------- labels (increment 3)
    post('labels_all', 'masked',
         {'patterns': [p(NDM_F, ident='NDM-F'), p(ABSENT_20, ident='absent')],
          'output': {'labels': 'all'}}, 200,
         'output.labels "all" within the threshold: every context with its labels (column, '
         'support kmer) and their placed occurrences (seq_id, record, strand, 1-based '
         'nt_coords, nt_length; the record mapping of the .seqs first, the offset after), '
         'by_label (contexts desc, column asc) with the deduplicated occurrences, counts.labels '
         'and counts.occurrences exact, placement record, annotation budgeted; the absent '
         'primer: complete, exact zeros',
         entries(expect_all(exact(24), complete, labelled, counted('labels', 'exact', 9),
                            counted('occurrences', 'exact', 42), field('placement', 'record')),
                 expect_all(exact(0), complete, counted('labels', 'exact', 0),
                            counted('occurrences', 'exact', 0)))),
    post('labels_all_withheld', 'masked',
         {'patterns': [p(NDM_F)], 'max_contexts': 10, 'output': {'labels': 'all'}}, 200,
         'output.labels "all" above the threshold: the count is not admitted, no annotation '
         'row is read (work.annotation_rows 0), withheld count_above_threshold, labels and '
         'occurrences unknown, by_label null',
         entries(expect_all(exact(24), withheld('count_above_threshold'),
                            counted('labels', 'unknown'), field('by_label', None)))),
    post('labels_all_truncated', 'masked',
         {'patterns': [p(NDM_F)], 'max_labels_per_anchor': 2, 'output': {'labels': 'all'}},
         200,
         'max_labels_per_anchor 2 below the rows\' 9 labels: all_or_count withholds '
         '(anchor_labels_truncated) and lists each truncated anchor (kmer, row, cap, total) '
         'so that the cap to ask for is known',
         entries(expect_all(exact(24), withheld('anchor_labels_truncated'),
                            lambda e: check(len(e['anchors_truncated']) == 24,
                                            e['anchors_truncated'])))),
    post('labels_all_truncated_partial', 'masked',
         {'patterns': [p(NDM_F)], 'mode': 'partial', 'max_labels_per_anchor': 2,
          'output': {'labels': 'all'}}, 200,
         'the same cap in mode partial: every context returned, its labels_status '
         '"truncated" with labels_total 9 and the first 2 labels (ascending column), the '
         'counts at_least, retrieval_complete false',
         entries(expect_all(counted('labels', 'at_least'), field('retrieval_complete', False),
                            lambda e: check(all(r['labels_status'] == 'truncated'
                                                for r in e['results']), e['results'])))),
    post('labels_all_partial', 'masked',
         {'patterns': [p(NDM_F)], 'mode': 'partial', 'max_contexts': 3, 'max_labels': 2,
          'max_occurrences_per_label': 1, 'output': {'labels': 'all'}}, 200,
         'mode partial with the label caps: the first 3 contexts (cut max_contexts), the first '
         '2 labels in label order (labels_cut), each label\'s first occurrence of its union '
         '(occurrences_cut); the counts over the returned contexts, at_least',
         entries(expect_all(cut('max_contexts'), field('returned', 3),
                            lambda e: check(e['labels_cut']['reason'] == 'max_labels'
                                            and e['occurrences_cut'] is not None, e)))),
    post('labels_all_work_budget', 'masked',
         {'patterns': [p(NDM_F), p(NDM_R)], 'max_annotation_work': 1,
          'output': {'labels': 'all'}}, 200,
         'max_annotation_work 1: the first read passes the budget, the reads stop (stop '
         '{label_discovery, max_annotation_work}), all_or_count withholds (annotation_budget); '
         'the budget is the request\'s: the next pattern reads nothing',
         entries(expect_all(exact(24), withheld('annotation_budget'),
                            field('stop', {'phase': 'label_discovery',
                                           'reason': 'max_annotation_work'})),
                 expect_all(withheld('annotation_budget'),
                            lambda e: check(e['work']['annotation_rows'] == 0, e['work'])))),
    post('labels_all_work_budget_partial', 'masked',
         {'patterns': [p(NDM_F)], 'mode': 'partial', 'max_annotation_work': 1,
          'output': {'labels': 'all'}}, 200,
         'the same stop in mode partial: every context returned, the row read first with its '
         'labels, the others labels_status "not_read" (labels null); counts at_least, '
         'retrieval_complete false',
         entries(expect_all(field('returned', 24), counted('labels', 'at_least'),
                            field('stop', {'phase': 'label_discovery',
                                           'reason': 'max_annotation_work'}),
                            lambda e: check(any(r['labels_status'] == 'not_read'
                                                for r in e['results']), e['results'])))),
    post('labels_all_global', 'masked_no_map',
         {'patterns': [p(NDM_F)], 'output': {'labels': 'all'}}, 200,
         'coordinates without the record mapping (--no-coord-mapping): placement global, each '
         'label\'s occurrence_list holds (kmer_coord, offset, strand), nothing is placed in a '
         'record (counts.occurrences unknown), note record_bounds_unknown',
         entries(expect_all(exact(24), complete, field('placement', 'global'),
                            field('notes', ['record_bounds_unknown']),
                            counted('occurrences', 'unknown')))),
    # GPT review 1 of the SPEC (2026-10-07), #11: the cases no fixture showed
    post('labels_all_partial_exact_cut', 'masked',
         {'patterns': [p(NDM_F)], 'mode': 'partial', 'max_labels': 2,
          'max_occurrences_per_label': 1, 'output': {'labels': 'all'}}, 200,
         'exact counts beside cut lists: mode partial returns every context (24, no cut) and '
         'reads every row, so counts.labels (9) and counts.occurrences (42) stay exact while '
         'max_labels 2 cuts by_label (labels_cut max_labels) and max_occurrences_per_label 1 '
         'each label\'s occurrence list (occurrences_cut); retrieval_complete false',
         entries(expect_all(exact(24), field('returned', 24), field('cut', None),
                            counted('labels', 'exact', 9), counted('occurrences', 'exact', 42),
                            field('retrieval_complete', False),
                            lambda e: check(e['labels_cut']['reason'] == 'max_labels'
                                            and e['occurrences_cut'] is not None, e)))),
    post('labels_all_mixed_slots', 'masked',
         {'patterns': [p(NDM_F, ident='ok'), p('ACGUACGTACGTACGTACGT', ident='bad'),
                       p(NDM_F11, ident='floor'), p(ABSENT_20, ident='absent')],
          'output': {'labels': 'all'}}, 200,
         'successful and refused slots in one retrieval request (all_or_count, labels "all"): '
         'the error slots (bad_alphabet, information_below_floor) carry only id, kind and '
         'error, spend nothing, and the slots around them are answered in full, labelled',
         entries(expect_all(exact(24), complete, labelled, counted('labels', 'exact', 9)),
                 slot_error('bad_alphabet'), slot_error('information_below_floor'),
                 expect_all(exact(0), complete, counted('labels', 'exact', 0)))),
    # (the two below rest on the memory model of SPEC
    # §14.4; their expects fail loudly if the account stops elsewhere. A refused row needs the
    # account nearly spent: on the mini only a request whose earlier patterns hold most of it)
    post('labels_all_output_budget', 'masked',
         {'patterns': [p('GCGGCGGCGGCG', 'dna', 'repeat')], 'max_memory_mb': 1,
          'max_contexts': 10000, 'output': {'labels': 'all'}}, 200,
         'the memory account (max_memory_mb 1, the smallest) holds the 1,828 contexts of a '
         'GCG repeat but not the labels and placed occurrences built for them: stop {output, '
         'max_memory}, all_or_count withholds (output_budget), the contexts count exact',
         entries(expect_all(exact(1828), withheld('output_budget'),
                            field('stop', {'phase': 'output', 'reason': 'max_memory'})))),
    post('labels_all_rows_refused', 'masked',
         {'patterns': [p('GCGGCGGCGGCG'), p('CGCCAGCGCCAG'), p('CGCCAGCGCCAG'), p('GCGGCGGCGGCG')],
          'mode': 'partial', 'max_memory_mb': 1, 'max_contexts': 10000,
          'output': {'labels': 'all'}}, 200,
         'one memory account over four patterns (max_memory_mb 1): the earlier patterns hold most '
         'of it, so a row of the last cannot be held alone and is refused: listed in rows_refused '
         '(phase, reason max_memory, needed_bytes), its contexts labels_status "refused", partial '
         'answers on; the labels count at_least',
         entries(lambda e: None, lambda e: None, lambda e: None,
                 expect_all(lambda e: check(len(e['rows_refused']) > 0, e['rows_refused']),
                            lambda e: check(any(r['labels_status'] == 'refused'
                                                for r in e['results']), e['results']),
                            counted('labels', 'at_least')))),
    post('labels_all_count', 'masked',
         {'patterns': [p(NDM_F)], 'mode': 'count', 'output': {'labels': 'all'}}, 200,
         'mode count with output.labels "all": counted, no annotation read, output null, the '
         'note annotation_not_read says the projection had no effect',
         entries(expect_all(exact(24), field('notes', ['annotation_not_read'])))),
    post('annotation_unbudgeted', 'primary',
         {'patterns': [p(NDM_F)], 'output': {'labels': 'all'}}, 400,
         '400 annotation_unbudgeted: output.labels "all" on an annotation without the '
         'budget-aware decode (a column annotation), unless the request allows it',
         refused('annotation_unbudgeted')),
    post('labels_all_unbudgeted', 'primary',
         {'patterns': [p(NDM_F)], 'output': {'labels': 'all'},
          'allow_unbudgeted_annotation': True}, 200,
         'the same with allow_unbudgeted_annotation: labels read without a memory bound on the '
         'reads (annotation unbudgeted, note annotation_unbudgeted); a PRIMARY index places '
         'nothing (placement none_canonical)',
         entries(expect_all(exact(24), complete, field('annotation', 'unbudgeted'),
                            field('placement', 'none_canonical')))),

    # ---------------------------------------------------------------- whole-request refusals
    post('unknown_field', 'masked',
         {'patterns': [p(NDM_F)], 'mode': 'count', 'bogus': 1}, 400,
         '400 invalid_request: a field nothing reads is refused, never ignored',
         refused('invalid_request')),
    post('too_many_patterns', 'masked',
         {'patterns': [p(NDM_F)] * 17, 'mode': 'count'}, 400,
         '400 invalid_request: more patterns than the server\'s max_patterns (16): a list is '
         'refused, never cut',
         refused('invalid_request')),
    post('later_increment_labels', 'masked',
         {'patterns': [p(NDM_F)], 'output': {'labels': 'predicate_only'}}, 400,
         '400 later_increment: output.labels "predicate_only" (the labels a predicate names) '
         'is not served yet (capabilities projections_later_increment; "all" is served since '
         'increment 3)',
         refused('later_increment')),
    post('later_increment_graphs', 'masked',
         {'patterns': [p(NDM_F)], 'mode': 'count', 'graphs': [INDEX_NAME]}, 400,
         '400 later_increment: the request field graphs (multi-graph selection) is refused by '
         'name, whatever its value',
         refused('later_increment')),
    post('resident_only', 'masked',
         {'patterns': [p(NDM_F)], 'mode': 'count', 'in_ram': False}, 400,
         '400 resident_only: in_ram, whatever its value (the route never loads an index)',
         refused('resident_only')),
    post('mask_required', 'unmasked',
         {'patterns': [p(NDM_F)], 'mode': 'count'}, 400,
         '400 mask_required: the graph has no dummy-edge mask, whatever the request asks',
         refused('mask_required')),
    post('multi_graph', 'multi',
         {'patterns': [p(NDM_F)], 'mode': 'count'}, 400,
         '400 later_increment on a multi-graph server, whatever the request',
         refused('later_increment')),
    post('representation_unsupported', 'hash',
         {'patterns': [p(NDM_F)], 'mode': 'count'}, 400,
         '400 representation_unsupported: a graph the engine does not recognise (a hash graph), '
         'whatever the request asks',
         refused('representation_unsupported')),

    # ---------------------------------------------------------------- hand-made
    post('deadline_503', 'masked',
         {'patterns': [p(NDM_F)], 'mode': 'partial', 'time_budget_ms': TINY_BUDGET_MS,
          'output': {'labels': 'none'}}, 503,
         'HAND-MADE: 503 deadline, the answer could not be written by time_budget_ms '
         '(nothing partial is sent)',
         hand_made={'error': 'pattern: the answer could not be written within time_budget_ms '
                             '(%d ms, the finalisation reserve of %d ms included): nothing '
                             'partial is sent' % (TINY_BUDGET_MS, FINALIZE_MS),
                    'code': 'deadline'}),
    post('initializing_503', 'masked',
         {'patterns': [p(NDM_F)], 'mode': 'count'}, 503,
         'HAND-MADE: 503 while the index loads (every route; no code; Retry-After: 60)',
         hand_made={'error': 'Server is currently initializing, please come back later.'},
         headers={'Retry-After': '60'}),
]


# ----------------------------------------------------------------------- servers

SERVERS = {
    'masked': 'server_query -i {work}/masked/graph_k31.dbg -a {work}/masked/' + MINI_ANNO
              + ' --index-name ' + INDEX_NAME + ' --index-release ' + INDEX_RELEASE
              + '  (a copy of the mini graph with the .edgemask of transform --mask-dummy; '
                'annotation, .seqs and sidecars linked from the mini)',
    'unmasked': 'server_query -i {mini}/graph_k31.dbg -a {mini}/' + MINI_ANNO
                + ' --index-name ' + INDEX_NAME + ' --index-release ' + INDEX_RELEASE
                + '  (the mini index as built: no .edgemask)',
    'built_at_load': 'server_query -i {mini}/graph_k31.dbg -a {mini}/' + MINI_ANNO
                     + ' --pattern-build-mask --index-name ' + INDEX_NAME + ' --index-release '
                     + INDEX_RELEASE + '  (the mini index as built, its mask built in memory at '
                       'start-up)',
    'multi': 'server_query {work}/graphs.csv --index-release ' + INDEX_RELEASE
             + '  (one line: ' + MULTI_GRAPH_NAME + ',{work}/masked/graph_k31.dbg,'
               '{work}/masked/' + MINI_ANNO + ')',
    'primary': 'server_query -i {work}/primary/graph.dbg -a {work}/primary/anno.column.annodbg'
               ' --index-name ' + PRIMARY_NAME + ' --index-release ' + INDEX_RELEASE
               + '  (a PRIMARY graph of ' + ' and '.join(PRIMARY_RECORDS)
               + ' at k = 31 with --mask-dummy, column annotation by file name)',
    'masked_no_map': 'server_query -i {work}/masked/graph_k31.dbg -a {work}/masked/' + MINI_ANNO
                     + ' --no-coord-mapping --index-name ' + INDEX_NAME + ' --index-release '
                     + INDEX_RELEASE + '  (the masked copy, its .seqs not loaded)',
    'hash': 'server_query -i {work}/hash/graph.orhashdbg -a {work}/hash/anno.column.annodbg'
            ' --index-release ' + INDEX_RELEASE + '  (a hash graph of ' + HASH_RECORD
            + ' at k = 31, column annotation by file name)',
}


def free_port():
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(('127.0.0.1', 0))
        return s.getsockname()[1]


def http(port, method, path, body=None):
    data = None if body is None else json.dumps(body).encode()
    req = urllib.request.Request(f'http://127.0.0.1:{port}{path}', data=data, method=method)
    try:
        with urllib.request.urlopen(req, timeout=300) as r:
            return r.status, r.read()
    except urllib.error.HTTPError as e:
        return e.code, e.read()


class Server:
    def __init__(self, binary, args, log_path):
        for _ in range(10):
            self.port = free_port()
            self.log = open(log_path, 'ab')
            self.process = subprocess.Popen(
                [binary, 'server_query'] + args
                + ['--port', str(self.port), '--address', '127.0.0.1', '-p', '2'],
                stdout=self.log, stderr=subprocess.STDOUT)
            if self._wait_ready(300):
                return
            self.stop()
        raise RuntimeError(f'could not start {binary} server_query {args} (log: {log_path})')

    def _wait_ready(self, timeout):
        deadline = time.time() + timeout
        while time.time() < deadline:
            if self.process.poll() is not None:
                return False
            try:
                if http(self.port, 'GET', '/stats')[0] == 200:
                    return True
            except (urllib.error.URLError, ConnectionError, OSError):
                pass
            time.sleep(0.2)
        return False

    def stop(self):
        if self.process.poll() is None:
            self.process.kill()
            self.process.wait()
        self.log.close()


def run(cmd, log):
    with open(log, 'ab') as f:
        f.write(('$ ' + ' '.join(cmd) + '\n').encode())
        f.flush()
        res = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
    if res.returncode:
        raise RuntimeError(f'failed ({res.returncode}): {" ".join(cmd)} (log: {log})')


def build_indexes(binary, mini, work):
    """The masked copy of the mini graph and the PRIMARY index (reused when present)."""
    log = os.path.join(work, 'build.log')
    fastas = sorted(glob.glob(os.path.join(mini, 'fasta', '*.fa')))
    if not fastas:
        raise RuntimeError(f'no FASTA under {mini}/fasta (scripts/traversal/build_mini_refseq.sh)')

    masked = os.path.join(work, 'masked')
    graph = os.path.join(masked, MINI_GRAPH)
    if not os.path.isfile(graph[:-len('.dbg')] + '.edgemask'):
        os.makedirs(masked, exist_ok=True)
        # the way a host gets its mask (DESIGN §4): transform writes only <graph>.edgemask beside
        # a copy of the graph, the mask `build --mask-dummy` writes (test_pattern.py and
        # test_transform_mask.py compare the two byte for byte); a copy, since the mini index
        # is only read
        shutil.copyfile(os.path.join(mini, MINI_GRAPH), graph)
        run([binary, 'transform', '--mask-dummy', '-p', '1', graph], log)
    if not filecmp.cmp(graph, os.path.join(mini, MINI_GRAPH), shallow=False):
        # node ids would differ: the mini's annotation would not describe this graph
        raise RuntimeError(f'{graph} differs from the mini index graph: rebuild the mini index')
    for f in [MINI_ANNO] + MINI_SIDECARS:
        link = os.path.join(masked, f)
        if not os.path.lexists(link):
            os.symlink(os.path.join(os.path.abspath(mini), f), link)

    primary = os.path.join(work, 'primary')
    if not os.path.isfile(os.path.join(primary, 'anno.column.annodbg')):
        os.makedirs(primary, exist_ok=True)
        for f in PRIMARY_RECORDS:
            shutil.copy(os.path.join(mini, 'fasta', f), os.path.join(primary, f))
        # relative paths (cwd primary/): the column labels are the file names
        def here(cmd):
            with open(log, 'ab') as f:
                f.write(('$ (cd primary) ' + ' '.join(cmd) + '\n').encode())
                f.flush()
                res = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, cwd=primary)
            if res.returncode:
                raise RuntimeError(f'failed: {" ".join(cmd)} (log: {log})')
        here([binary, 'build', '-p', '1', '--mode', 'canonical', '--graph', 'succinct',
              '-k', str(MINI_K), '--in-ram', '-o', 'canonical'] + PRIMARY_RECORDS)
        here([binary, 'transform', '-p', '1', '--to-fasta', '--primary-kmers', '-o', 'primary',
              'canonical.dbg'])
        here([binary, 'build', '-p', '1', '--mode', 'primary', '--graph', 'succinct',
              '-k', str(MINI_K), '--mask-dummy', '--in-ram', '-o', 'graph',
              'primary.fasta.gz'])
        here([binary, 'annotate', '-p', '1', '-i', 'graph.dbg', '--anno-filename', '-o', 'anno']
             + PRIMARY_RECORDS)

    hashed = os.path.join(work, 'hash')
    if not os.path.isfile(os.path.join(hashed, 'anno.column.annodbg')):
        os.makedirs(hashed, exist_ok=True)
        shutil.copy(os.path.join(mini, 'fasta', HASH_RECORD), os.path.join(hashed, HASH_RECORD))
        # a graph the engine does not recognise: neither a DBGSuccinct nor a wrapped PRIMARY one
        for cmd in ([binary, 'build', '-p', '1', '--graph', 'hash', '-k', str(MINI_K),
                     '-o', 'graph', HASH_RECORD],
                    [binary, 'annotate', '-p', '1', '-i', 'graph.orhashdbg', '--anno-filename',
                     '-o', 'anno', HASH_RECORD]):
            with open(log, 'ab') as f:
                f.write(('$ (cd hash) ' + ' '.join(cmd) + '\n').encode())
                f.flush()
                res = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, cwd=hashed)
            if res.returncode:
                raise RuntimeError(f'failed: {" ".join(cmd)} (log: {log})')

    with open(os.path.join(work, 'graphs.csv'), 'w') as f:
        f.write(f'{MULTI_GRAPH_NAME},{graph},{os.path.join(masked, MINI_ANNO)}\n')


def server_args(name, mini, work):
    masked = os.path.join(work, 'masked')
    ident = ['--index-name', INDEX_NAME, '--index-release', INDEX_RELEASE]
    if name == 'masked':
        return ['-i', os.path.join(masked, MINI_GRAPH), '-a', os.path.join(masked, MINI_ANNO)] \
            + ident
    if name == 'unmasked':
        return ['-i', os.path.join(mini, MINI_GRAPH), '-a', os.path.join(mini, MINI_ANNO)] \
            + ident
    if name == 'built_at_load':
        return ['-i', os.path.join(mini, MINI_GRAPH), '-a', os.path.join(mini, MINI_ANNO),
                '--pattern-build-mask'] + ident
    if name == 'multi':
        # one name and one manifest describe one index: a multi-graph server takes neither
        return [os.path.join(work, 'graphs.csv'), '--index-release', INDEX_RELEASE]
    if name == 'masked_no_map':
        return ['-i', os.path.join(masked, MINI_GRAPH), '-a', os.path.join(masked, MINI_ANNO),
                '--no-coord-mapping'] + ident
    if name == 'hash':
        hashed = os.path.join(work, 'hash')
        return ['-i', os.path.join(hashed, 'graph.orhashdbg'),
                '-a', os.path.join(hashed, 'anno.column.annodbg'),
                '--index-release', INDEX_RELEASE]
    if name == 'primary':
        primary = os.path.join(work, 'primary')
        return ['-i', os.path.join(primary, 'graph.dbg'),
                '-a', os.path.join(primary, 'anno.column.annodbg'),
                '--index-name', PRIMARY_NAME, '--index-release', INDEX_RELEASE]
    raise ValueError(name)


# ----------------------------------------------------------------------- files

def dumps(value):
    return json.dumps(value, indent=1, ensure_ascii=False) + '\n'


def replace_paths(value, prefixes):
    """Every string under a work directory written as {work}/... (the mini as {mini}/...)."""
    if isinstance(value, dict):
        return {k: replace_paths(v, prefixes) for k, v in value.items()}
    if isinstance(value, list):
        return [replace_paths(v, prefixes) for v in value]
    if isinstance(value, str):
        for prefix, name in prefixes:
            value = value.replace(prefix, name)
    return value


BLANK = '<varies>'


def blank_count(c):
    if isinstance(c, dict) and 'relation' in c and 'unit' in c:
        out = {k: (blank_count(v) if isinstance(v, dict) else v) for k, v in c.items()
               if k not in ('lower', 'upper')}
        out['value'] = BLANK
        out['relation'] = BLANK
        return out
    if isinstance(c, dict):
        return {k: blank_count(v) for k, v in c.items()}
    return c


def blanked(answer):
    """The answer without what varies between runs (see the module docstring)."""
    a = copy.deepcopy(answer)

    def walk(x):
        if isinstance(x, dict):
            for k in x:
                if k == 'server_instance':
                    x[k] = BLANK
                else:
                    walk(x[k])
        elif isinstance(x, list):
            for v in x:
                walk(v)
    walk(a)
    if isinstance(a, dict) and isinstance(a.get('patterns'), list):
        if isinstance(a.get('timing'), dict):
            a['timing']['elapsed_ms'] = BLANK
        clocked = False
        for e in a['patterns']:
            if isinstance(e.get('timing'), dict):
                # elapsed_ms, and with output.labels "all" label_discovery_ms, placement_ms
                for k in e['timing']:
                    e['timing'][k] = BLANK
            if e.get('determinism') == 'time_limited' and not clocked:
                # where the clock stopped the search: how far it got. Only the first such
                # pattern: the later ones are stopped by the same budget before they start,
                # the same in every run (review of 2026-10-07, C2-02)
                clocked = True
                e['counts'] = blank_count(e['counts'])
                e['work'] = {k: BLANK for k in e['work']}
                if e.get('cut') == {'reason': 'time'} \
                        and e.get('stop') != {'phase': 'discovery', 'reason': 'time'}:
                    # the clock cut its release: how many contexts got out (a time stop in
                    # discovery releases nothing, returned 0, compared; review of 2026-10-07,
                    # C1-03)
                    e['returned'] = BLANK
                    e['results'] = BLANK
    return a


def read_text(path):
    if not os.path.isfile(path):
        return None
    with open(path, encoding='utf-8') as f:
        return f.read()


def readme(fixtures):
    lines = [
        '# POST /pattern fixtures (pattern_contract_version 1)',
        '',
        'The bodies of `docs/SPEC-pattern-search.md`, generated by',
        '`scripts/traversal/pattern_fixtures.py` from a server_query on copies of the mini index',
        '(`build/mini_refseq`): one directory per fixture with `request.json` (the body sent;',
        '`null` for a GET) and `answer.json` (the body answered). `index.json` states each',
        'fixture\'s server, method, path and HTTP status, and the servers\' command lines.',
        'Regenerate with `pattern_fixtures.py`, verify with `pattern_fixtures.py --check`;',
        '`api/python/tests/test_pattern_fixtures.py` validates every answer against the SPEC\'s',
        'field lists.',
        '',
        'Varies between runs (stored as answered, blanked by --check): `timing.elapsed_ms`',
        '(and every `timing` value of a pattern), `server_instance`, and the counts and work of',
        'the first `determinism: time_limited` pattern of an answer (and, when the clock cut its',
        'release, its `returned` and `results`; the later ones, stopped by the same budget, are',
        'compared as they are).',
        'Paths under the generator\'s work directory read `{work}/...`. Two fixtures are',
        'HAND-MADE (no server produces them on demand); their text is the code\'s. No fixture',
        'holds `primary_unwrapped` (server_query and the CLI always wrap a PRIMARY graph),',
        '`alphabet_unsupported` (a graph of another alphabet does not load in this build),',
        '`mask_invalid` (a graph extended after masking by an older build, or a stale mask) or',
        '`alphabet_untested` (a DNA5 build); test_pattern_fixtures.py names them in NO_FIXTURE.',
        '',
    ]
    for f in fixtures:
        lines.append(f'- `{f.name}` ({f.method} {f.path.split("?")[0]}, {f.status}, '
                     f'{f.server}): {f.shows}')
    return '\n'.join(lines) + '\n'


def index_json(fixtures):
    return {
        'pattern_contract_version': 1,
        'generator': 'scripts/traversal/pattern_fixtures.py',
        'spec': 'docs/SPEC-pattern-search.md',
        'servers': SERVERS,
        'fixtures': {f.name: {'server': f.server, 'method': f.method, 'path': f.path,
                              'status': f.status, 'hand_made': f.hand_made is not None,
                              'headers': f.headers, 'shows': f.shows}
                     for f in fixtures},
    }


# ----------------------------------------------------------------------- main

def generate(args, fixtures):
    """{name: (request, answer)} as the servers answer now."""
    out = {}
    for f in fixtures:
        if f.hand_made is not None:
            out[f.name] = (f.body, f.hand_made)
    by_server = {}
    for f in fixtures:
        if f.hand_made is None:
            by_server.setdefault(f.server, []).append(f)
    if not by_server:
        return out

    work = args.work or tempfile.mkdtemp(prefix='pattern_fixtures.')
    os.makedirs(work, exist_ok=True)
    work = os.path.realpath(work)
    mini = os.path.realpath(args.mini)
    prefixes = [(work, '{work}'), (mini, '{mini}')]
    try:
        build_indexes(args.metagraph, mini, work)
        for name, todo in by_server.items():
            server = Server(args.metagraph, server_args(name, mini, work),
                            os.path.join(work, f'server_{name}.log'))
            try:
                for f in todo:
                    status, raw = http(server.port, f.method, f.path, f.body)
                    answer = replace_paths(json.loads(raw), prefixes)
                    if status != f.status:
                        raise AssertionError(f'{f.name}: HTTP {status}, expected {f.status}: '
                                             f'{raw[:500]!r}')
                    if f.expect:
                        try:
                            f.expect(answer)
                        except (AssertionError, KeyError, TypeError) as e:
                            raise AssertionError(f'{f.name} no longer shows what it is for '
                                                 f'({f.shows}): {e!r}') from e
                    out[f.name] = (f.body, answer)
            finally:
                server.stop()
    finally:
        if not args.work:
            shutil.rmtree(work, ignore_errors=True)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--metagraph', default=DEFAULT_BIN)
    parser.add_argument('--mini', default=DEFAULT_MINI)
    parser.add_argument('--dir', default=DEFAULT_DIR)
    parser.add_argument('--work', default=None)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--only', nargs='+', default=None)
    args = parser.parse_args()
    # the PRIMARY index is built with its directory as the working directory (its column
    # labels are file names), where a relative binary path would name nothing
    args.metagraph = os.path.abspath(args.metagraph)

    names = [f.name for f in FIXTURES]
    assert len(set(names)) == len(names), 'fixture names must be unique'
    fixtures = [f for f in FIXTURES if not args.only or f.name in args.only]
    unknown = set(args.only or []) - set(names)
    if unknown:
        parser.error(f'unknown fixtures: {sorted(unknown)}')

    bodies = generate(args, fixtures)

    differ = []
    for f in fixtures:
        request, answer = bodies[f.name]
        d = os.path.join(args.dir, f.name)
        req_path, ans_path = os.path.join(d, 'request.json'), os.path.join(d, 'answer.json')
        stored = None
        if os.path.isfile(ans_path):
            with open(ans_path, encoding='utf-8') as fh:
                stored_text = fh.read()
            stored = json.loads(stored_text)
            # the stored file must be in canonical form, and equal once blanked
            same_answer = (stored_text == dumps(stored)
                           and dumps(blanked(stored)) == dumps(blanked(answer)))
        else:
            same_answer = False
        same_request = read_text(req_path) == dumps(request)
        if args.check:
            if not same_request:
                differ.append(req_path)
            if not same_answer:
                differ.append(ans_path)
            continue
        os.makedirs(d, exist_ok=True)
        if not same_request:
            with open(req_path, 'w', encoding='utf-8') as fh:
                fh.write(dumps(request))
        if not same_answer:
            with open(ans_path, 'w', encoding='utf-8') as fh:
                fh.write(dumps(answer))
            print(f'wrote {ans_path}')

    if not args.only:
        for path, text in ((os.path.join(args.dir, 'README.md'), readme(FIXTURES)),
                           (os.path.join(args.dir, 'index.json'), dumps(index_json(FIXTURES)))):
            if read_text(path) == text:
                continue
            if args.check:
                differ.append(path)
            else:
                with open(path, 'w', encoding='utf-8') as fh:
                    fh.write(text)
        stale = sorted(n for n in (os.listdir(args.dir) if os.path.isdir(args.dir) else [])
                       if os.path.isdir(os.path.join(args.dir, n)) and n not in names)
        if stale:
            if args.check:
                differ.extend(os.path.join(args.dir, n) for n in stale)
            else:
                for n in stale:
                    shutil.rmtree(os.path.join(args.dir, n))
                    print(f'removed stale fixture {n}')

    if args.check:
        if differ:
            print('differs from what the server answers now:\n  ' + '\n  '.join(differ))
            return 1
        print(f'OK: {len(fixtures)} fixtures match')
    return 0


if __name__ == '__main__':
    sys.exit(main())
