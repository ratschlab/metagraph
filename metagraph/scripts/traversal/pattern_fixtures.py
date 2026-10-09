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
  unmasked   the mini index as built (no .edgemask), served all the same (capabilities counting
             upper_bound): a count the search could not resolve is the bounds [lower, U], U the
             graph's candidate entries (its source dummies among them), with the additive
             estimate U x f (f the dummy fraction sampled at load, index.dummy_fraction); the
             lists stay exact. A pattern with at most --pattern-max-checked-entries (default
             50) unchecked candidates has each tested at query time, its counts then exact
  unmasked_unchecked
             the same with --pattern-max-checked-entries 0: no candidate is tested, so every
             count with unchecked candidates is the bounds [lower, U] -- what any pattern above
             the limit answers on the default server, shown on numbers small enough to follow
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
  masked_small_predicate_cap
             the masked copy served with --pattern-max-predicate-labels 4: a predicate of more
             names is refused (predicate_too_large)
  hash       a hash graph (build --graph hash) of one mini record (1296536.fa) with a column
             annotation: a graph the engine does not recognise (representation_unsupported)

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
same in every run, and compared as they are. Paths
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
# steps on the mini (about 1.5 s): stopped by any budget of a few milliseconds. (Not 12 Ns then
# NDM-F's first 12 bases: the engine skips a pattern's leading and trailing N runs, so that one
# takes 125 steps)
HEAVY = NDM_F[:2] + 'N' * 19 + NDM_F[3:13]
# a 10-mer whose suffix-scope count needs a mask scan (the prefix of a record's first k-mer):
# discovery takes 20 range steps, its scan one more, so max_steps 20 stops in the scan with
# bounds [lower, upper] (a body that carries bounds and mask_scan)
SCAN_10 = 'ATGCCGGTGA'
SCAN_STEPS = 20
# the record the hash graph is built from (a graph the engine does not recognise)
HASH_RECORD = '1296536.fa'
# graphs without the dummy-edge mask (DESIGN-pattern-search.md §4.4): the first 16 bases of an
# E. coli record (562.fa, NZ_CP021206.1) whose first k-mer no k-mer enters, so that the graph's
# source dummies ($^j and its first k - j bases, 1 <= j <= 15) hold it at offsets 1 to 15: 17
# contexts on the masked graph (1 on +, 16 on -), an upper bound of 32 without the mask (15
# source dummies more), its estimate 32 -- the estimate is not a bound
START16 = 'GATGCCGGTGAACAAC'
START16_EXACT = 17
START16_UPPER = 32
# a few unchecked candidates are tested one by one (max_checked_entries,
# DESIGN-pattern-search.md §4.4): the first 14 bases of the same record, [4, 72] without the
# check: 68 unchecked candidates, more than the default limit of 50, so it stays bounds on the
# default server too
START14 = START16[:14]
START14_LOWER = 4
START14_UPPER = 72
# paths (long_search "paths"): a 51-mer of two copies of a repeated 31-mer of the
# E. coli records (562): 10 bases before one copy, the 31-mer, 10 bases after another. Every
# k-mer of it is in the index (each window lies in one copy), so it is a path of the graph, but
# no record holds it whole: its labels are carried (label_intersection), none record_verified
CHIMERA = 'CGGCCTCCAGAGCACTTTGTCGTTTTTGGACGGAAAATCCCTAGAACCCCT'
# peptides: the first residues of the NDM-1 protein (blaNDM-1's first codons,
# ATG GAA TTG CCC AAT ATT ATG CAC CCG GTC GCG AAG CTG AGC): 10 residues are 30 bases, within
# one k-mer; 14 residues are 42 bases, a pattern longer than k
NDM_PEP = 'MELPNIMHPV'
NDM_PEP_LONG = 'MELPNIMHPVAKLS'
# '*' is a stop codon of the request's genetic code. NDM-1's last 9 residues
# and its stop (... ACG GCC CGC ATG GCC GAC AAG CTG CGC TGA): 30 bases, within one k-mer
NDM_PEP_STOP = 'TARMADKLR*'
# the same with X (any residue, never a stop) in place of the stop: no instance on the mini
NDM_PEP_STOP_X = 'TARMADKLRX'

# predicates (SPEC §19): a GCG repeat with 1,828 contexts on the mini, 1,760 of
# them carrying 287 (Pseudomonas aeruginosa) only; the nine columns of the mini (taxids as
# strings); a selective predicate (an E. coli column, not P. aeruginosa's) that keeps 60 of them;
# a typo; and the server variant's cap on a predicate's names (masked_small_predicate_cap)
GCG12 = 'GCGGCGGCGGCG'
MINI_COLUMNS = ['1296536', '158836', '287', '470', '546', '562', '573', '615', '72407']
SELECTIVE = {'and': [{'any': ['562']}, {'none': ['287']}]}
SELECTIVE_COUNT = 60
TYPO = '5622'
SMALL_PREDICATE_CAP = 4

# the server's defaults (--pattern-* flags, DESIGN-pattern-search.md §5.3), for the hand-made
# 503 and the expectations
FINALIZE_MS = 250
TINY_BUDGET_MS = FINALIZE_MS + 1
# a step stop whose release then meets the work time: ACG in
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


def paths_count(relation, value=None, extension=None):
    """counts.paths (long_search "paths") with |relation|, |value| and what the
    extension did, when given."""
    def run(entry):
        c = entry['counts']['paths']
        check(c['relation'] == relation and (value is None or c['value'] == value)
              and (extension is None or c['extension'] == extension), c)
    return run


def label_supports(*supports):
    """Every label of every path result has a support of |supports|, and each is shown."""
    def run(entry):
        seen = {label['support'] for r in entry['results'] for label in r['labels']}
        check(seen == set(supports), seen)
    return run


def has_note(name):
    def run(entry):
        check(name in entry['notes'], entry['notes'])
    return run


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


def selection(pass_, tested=None, selected=None, access='rows'):
    """SPEC §19.7: selection.pass, counts.tested and counts.selected (each a
    (relation, value) pair, or for selected bounds (lower, upper)), and the access."""
    def run(entry):
        check(entry['selection']['pass'] == pass_ and entry['selection']['access'] == access,
              entry['selection'])
        check(entry['absence_filter'] == 'predicate', entry.get('absence_filter'))
        c = entry['counts']
        if tested is not None:
            check((c['tested']['relation'], c['tested']['value']) == tested, c['tested'])
        if selected is not None:
            got = (c['selected']['lower'], c['selected']['upper']) \
                if c['selected']['relation'] == 'bounds' \
                else (c['selected']['relation'], c['selected']['value'])
            check(got == selected, c['selected'])
    return run


def selection_labels_are(*names):
    """Every result's selection_labels is one of the lists |names|."""
    def run(entry):
        lists = {tuple(r['selection_labels']) for r in entry['results']}
        check(lists <= {tuple(n) for n in names} and entry['results'], lists)
    return run


def selection_strands_are(*strands):
    """Every result's selection_strands (per label the orientation whose row carries it) is one
    of the lists |strands|, one per selection label."""
    def run(entry):
        lists = {tuple(r['selection_strands']) for r in entry['results']}
        check(lists <= {tuple(n) for n in strands} and entry['results'], lists)
        check(all(len(r['selection_strands']) == len(r['selection_labels'])
                  for r in entry['results']), entry['results'])
    return run


def predicate_block(normal_form=None, unknown=None, vacuous=None, strands=None):
    """The answer's predicate block (SPEC §19.10)."""
    def run(answer):
        b = answer['predicate']
        if normal_form is not None:
            check(b['normal_form'] == normal_form, b)
        if unknown is not None:
            check(b['unknown_labels'] == unknown, b)
        if vacuous is not None:
            check(b['vacuous'] is vacuous, b)
        if strands is not None:
            check(b['strands'] == strands, b)
    return run


def caps_block(available, reason=None, mask=None, counting=None):
    def run(doc):
        b = doc['pattern']
        check(b['pattern_contract_version'] == 1, b)
        check(b['available'] is available and b['unavailable_reason'] == reason, b)
        if mask is not None:
            check(b['mask'] == mask, b)
        if counting is not None:
            # exact with the mask, upper_bound (and the dummy fraction the
            # estimates rest on) without it
            check(b['counting'] == counting, b)
            check((b['dummy_fraction'] is not None) is (counting == 'upper_bound'), b)
    return run


def estimated(lower, upper, estimate):
    """A graph without its mask: counts.contexts (or anchors) is the bounds
    [lower, upper] with this estimate, and the entry says so (estimate_sampled_dummy_fraction)."""
    def run(entry):
        c = counts_of(entry)
        check((c['relation'], c['lower'], c['upper'], c.get('estimate'))
              == ('bounds', lower, upper, estimate), c)
        check('estimate_sampled_dummy_fraction' in entry['notes'], entry['notes'])
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
        'features and routes, the block available (basic, mask file, counting exact, placement '
        'record)',
        expect_all(features(True), caps_block(True, mask='file', counting='exact'))),
    get('traverse_capabilities', 'masked', '/traverse/capabilities',
        'GET /traverse/capabilities (the document the service probe reads) on the same server: '
        'the same `pattern` block',
        caps_block(True, mask='file', counting='exact')),
    get('capabilities_mask_absent', 'unmasked', '/capabilities',
        'the mini index as built, without its .edgemask (owner decision #16): the feature and '
        'route listed, the block available, mask absent, counting upper_bound with the '
        'dummy_fraction sampled at load (value, 95% interval, samples, source sampled); it '
        'said available false, mask_required before',
        expect_all(features(True), caps_block(True, mask='absent', counting='upper_bound'))),
    get('traverse_capabilities_mask_absent', 'unmasked', '/traverse/capabilities',
        'the same unmasked server on the probe route',
        caps_block(True, mask='absent', counting='upper_bound')),
    get('capabilities_built_at_load', 'built_at_load', '/capabilities',
        'the mini index as built, served with --pattern-build-mask: the block available, mask '
        'built_at_load (the mask built in memory at start-up), counting exact, otherwise as '
        'with the file',
        expect_all(features(True), caps_block(True, mask='built_at_load', counting='exact'))),
    get('traverse_capabilities_built_at_load', 'built_at_load', '/traverse/capabilities',
        'the same server on the probe route',
        caps_block(True, mask='built_at_load', counting='exact')),
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

    # ---------------------------------------------------------------- labels
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
    # more cases of output.labels "all"
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

    # ---------------------------------------------------------------- paths
    post('paths', 'masked',
         {'patterns': [p(NDM_40, ident='NDM-40'), p(ABSENT_40, ident='absent')],
          'long_search': 'paths', 'output': {'labels': 'none'}}, 200,
         'long_search "paths" (increment 4, opt-in): a pattern longer than k extended into '
         'paths: counts.paths exact with its per-strand split, candidates_examined and extension '
         'completed, work.extension_edges, timing.extension_ms; each path result has sequence '
         '(the L bases), anchor_kmer, instance, offset 0, strand, nodes and rows (never kmer); '
         'limits echo long_search and max_paths; no anchor: no path, the empty answer complete',
         entries(expect_all(paths_count('exact', 2, 'completed'), complete, field('returned', 2)),
                 expect_all(paths_count('exact', 0, 'no_anchors'), complete))),
    post('paths_count', 'masked',
         {'patterns': [p(NDM_40)], 'mode': 'count', 'long_search': 'paths'}, 200,
         'mode count with long_search "paths": the paths counted (exact, by strand), nothing '
         'released',
         entries(paths_count('exact', 2, 'completed'))),
    post('paths_labels', 'masked',
         {'patterns': [p(NDM_40, ident='NDM-40'), p(CHIMERA, ident='chimera')],
          'long_search': 'paths', 'output': {'labels': 'all'}}, 200,
         'the labels of paths (owner decision #14): each label with its support -- '
         'record_verified where one record holds the whole path (blaNDM-1\'s first 40 bases: '
         '9 columns, their placed occurrences of the whole path), label_intersection where '
         'every k-mer of the path carries the label but no record holds it whole (a 51-mer '
         'joining two copies of a repeated 31-mer of the E. coli records: occurrences exact '
         '0); by_label with paths and paths_record_verified, counts.labels.by_support',
         entries(expect_all(paths_count('exact', 2), complete, label_supports('record_verified'),
                            counted('labels', 'exact', 9), counted('occurrences', 'exact', 42),
                            field('placement', 'record')),
                 expect_all(paths_count('exact', 2), complete,
                            label_supports('label_intersection'),
                            counted('occurrences', 'exact', 0)))),
    post('paths_require_support', 'masked',
         {'patterns': [p(NDM_40, ident='NDM-40'), p(CHIMERA, ident='chimera')],
          'long_search': 'paths', 'require_support': 'record_verified',
          'output': {'labels': 'all'}}, 200,
         'require_support "record_verified": only the labels one record verifies are listed; '
         'each path and the entry count the others in labels_excluded_unverified (the '
         'chimera\'s two columns), limits.require_support echoed',
         entries(expect_all(paths_count('exact', 2), complete, label_supports('record_verified'),
                            counted('labels', 'exact', 9)),
                 expect_all(paths_count('exact', 2), complete,
                            lambda e: check(all(r['labels'] == [] for r in e['results']), e),
                            lambda e: check(e['labels_excluded_unverified']['value'] == 2,
                                            e['labels_excluded_unverified'])))),
    post('paths_max_paths_partial', 'masked',
         {'patterns': [p(NDM_40)], 'mode': 'partial', 'long_search': 'paths', 'max_paths': 1,
          'output': {'labels': 'none'}}, 200,
         'partial, more paths (2) than max_paths (1): the first path in answer order, '
         'cut.reason max_paths, the count exact',
         entries(expect_all(paths_count('exact', 2), cut('max_paths'), field('returned', 1),
                            field('stop', None)))),
    post('paths_count_above_threshold', 'masked',
         {'patterns': [p(NDM_40)], 'long_search': 'paths', 'max_paths': 1,
          'output': {'labels': 'none'}}, 200,
         'all_or_count, the exact path count (2) above max_paths (1): withheld '
         'count_above_threshold',
         entries(expect_all(paths_count('exact', 2), withheld('count_above_threshold')))),
    post('paths_stop_at_max_paths', 'masked',
         {'patterns': [p(NDM_40)], 'long_search': 'paths', 'max_paths': 1,
          'stop_at_threshold': True, 'output': {'labels': 'none'}}, 200,
         'stop_at_threshold with max_paths 1: the extension stops once more than max_paths '
         'paths are complete: stop {extension, max_paths}, paths at_least, extension stopped, '
         'withheld threshold_crossed',
         entries(expect_all(paths_count('at_least', 2, 'stopped'), withheld('threshold_crossed'),
                            field('stop', {'phase': 'extension', 'reason': 'max_paths'})))),
    post('paths_anchors_above_threshold', 'masked',
         {'patterns': [p(NDM_40)], 'long_search': 'paths', 'max_anchors': 1,
          'output': {'labels': 'none'}}, 200,
         'the exact anchor count (2) above max_anchors (1): the extension is not admitted, '
         'paths unknown (extension not_admitted), withheld anchors_above_threshold',
         entries(expect_all(exact(2), paths_count('unknown', None, 'not_admitted'),
                            withheld('anchors_above_threshold')))),
    post('paths_max_steps', 'masked',
         {'patterns': [p(NDM_40, ident='NDM-40'), p(NDM_F, ident='NDM-F')], 'mode': 'partial',
          'long_search': 'paths', 'max_steps': 75, 'output': {'labels': 'none'}}, 200,
         'max_steps reached in the extension (discovery took 62 steps): stop {extension, '
         'max_steps}, paths at_least, the paths completed before the stop released (a prefix of '
         'the answer order), cut max_steps; the next pattern answers the sticky stop',
         entries(expect_all(paths_count('at_least', 1, 'stopped'), cut('max_steps'),
                            field('returned', 1),
                            field('stop', {'phase': 'extension', 'reason': 'max_steps'})),
                 expect_all(cut('max_steps'), field('returned', 0), unknown_counts,
                            field('stop', {'phase': 'discovery', 'reason': 'max_steps'})))),
    post('paths_global', 'masked_no_map',
         {'patterns': [p(NDM_40)], 'long_search': 'paths', 'output': {'labels': 'all'}}, 200,
         'paths with coordinates but no record mapping (placement global): nothing can be '
         'verified, every label label_intersection, its occurrence_list the chains '
         '(kmer_coord of the first k-mer, offset 0, strand), note record_bounds_unknown',
         entries(expect_all(paths_count('exact', 2), complete, label_supports('label_intersection'),
                            field('placement', 'global'), has_note('record_bounds_unknown')))),
    post('paths_primary', 'primary',
         {'patterns': [p(NDM_40)], 'long_search': 'paths', 'output': {'labels': 'all'},
          'allow_unbudgeted_annotation': True}, 200,
         'paths on a PRIMARY index: orientation instead of strand, wrapper node ids; no '
         'coordinates read (placement none_canonical): every label label_intersection, note '
         'label_intersection_only',
         entries(expect_all(paths_count('exact', 2), complete, label_supports('label_intersection'),
                            field('placement', 'none_canonical'),
                            has_note('label_intersection_only')))),
    post('support_unavailable', 'masked_no_map',
         {'patterns': [p(NDM_40)], 'long_search': 'paths', 'require_support': 'record_verified',
          'output': {'labels': 'all'}}, 400,
         '400 support_unavailable: require_support "record_verified" on an index that cannot '
         'verify (no record mapping: its best support is label_intersection)',
         refused('support_unavailable')),

    # ---------------------------------------------------------------- peptides
    post('peptide', 'masked',
         {'patterns': [p(NDM_PEP, 'protein', 'NDM-1 1-10'), p('MELPNJMHPV', 'protein', 'J'),
                       p(NDM_PEP_LONG, 'protein', 'NDM-1 1-14')],
          'output': {'labels': 'all'}}, 200,
         'protein patterns (increment 5) in the default genetic code (1): kind protein, length '
         'in bases (3 per residue), residues, genetic_code; the first 10 residues of NDM-1 (30 '
         'bases, one k-mer) in their codons: the contexts, labels and placed occurrences of '
         'blaNDM-1 (9 columns, 42); J (I or L) admits the same; 14 residues (42 bases, longer '
         'than k) without long_search: anchors only, withheld paths_later_increment',
         entries(expect_all(exact(4), complete, labelled, counted('labels', 'exact', 9),
                            counted('occurrences', 'exact', 42), field('residues', 10),
                            field('length', 30)),
                 expect_all(exact(4), complete, counted('occurrences', 'exact', 42)),
                 expect_all(exact(2), withheld('paths_later_increment'),
                            field('residues', 14)))),
    post('peptide_count', 'masked',
         {'patterns': [p(NDM_PEP, 'protein')], 'mode': 'count', 'genetic_code': 11}, 200,
         'a peptide counted in the bacterial genetic code (genetic_code 11, stated in the entry)',
         entries(expect_all(exact(4), field('genetic_code', 11)))),
    post('peptide_paths', 'masked',
         {'patterns': [p(NDM_PEP_LONG, 'protein', 'NDM-1 1-14')], 'long_search': 'paths',
          'output': {'labels': 'all'}}, 200,
         'a peptide longer than k with long_search "paths": its anchors extended through the '
         'codon automaton into paths of 42 bases (one per strand), each with its labels, '
         'record_verified',
         entries(expect_all(paths_count('exact', 2, 'completed'), complete,
                            label_supports('record_verified'), counted('labels', 'exact', 9),
                            counted('occurrences', 'exact', 42)))),
    post('peptide_bad_residue', 'masked',
         {'patterns': [p('MELPUIMHPV', 'protein', 'U'), p('MELPNIMHPV*', 'protein', 'stop'),
                       p(NDM_PEP, 'protein', 'ok')], 'mode': 'count'}, 200,
         'error slot bad_alphabet of a protein pattern (U, selenocysteine, is not served); the '
         'stop * is a residue since owner decision #19 (a stop codon of the genetic code): '
         'MELPNIMHPV* (33 bases, longer than k) is answered, its anchors counted (it was '
         'refused stop_unsupported by 4596bb3b, a code now retired); the last pattern is answered',
         entries(slot_error('bad_alphabet'), expect_all(field('scope', 'long'),
                                                        field('residues', 11)),
                 exact(4))),
    post('peptide_stop', 'masked',
         {'patterns': [p(NDM_PEP_STOP, 'protein', 'NDM-1 end'),
                       p(NDM_PEP_STOP_X, 'protein', 'X')], 'output': {'labels': 'none'}}, 200,
         'the stop * (owner decision #19) in the standard code (1): a stop codon (TAA, TAG, TGA) '
         'at that position; NDM-1\'s last 9 residues and its stop, 30 bases: the contexts end '
         'in its stop codon TGA; X (any residue) never matches a stop: the same peptide with X '
         'for * has no context',
         entries(expect_all(exact(4), complete,
                            lambda e: check(all(r['instance'].endswith('TGA') if r['strand'] == '+'
                                                else r['instance'].startswith('TCA')
                                                for r in e['results']), e['results'])),
                 expect_all(exact(0), complete))),
    post('peptide_no_stop_codon', 'masked',
         {'patterns': [p(NDM_PEP_STOP, 'protein', 'NDM-1 end'),
                       p(NDM_PEP_LONG + '*', 'protein', 'NDM-1 1-14 and a stop')],
          'mode': 'count', 'genetic_code': 27, 'long_search': 'paths'}, 200,
         'a table without an unconditional stop codon (27; also 28 and 31: their context stops '
         'match as their residue, owner decision #21): * matches nothing, and the answer says '
         'so (note no_stop_codon): 30 bases, exact 0 without a search (work 0); 45 bases, '
         'longer than k: its anchor window holds no *, its anchors counted, its paths exact 0',
         entries(expect_all(exact(0), has_note('no_stop_codon'),
                            lambda e: check(e['work']['steps'] == 0, e['work'])),
                 expect_all(has_note('no_stop_codon'), paths_count('exact', 0, 'completed')))),
    post('peptide_no_stop_codon_after_stop', 'masked',
         {'patterns': [p(HEAVY, 'iupac', 'heavy'), p(NDM_PEP_STOP, 'protein', 'NDM-1 end')],
          'mode': 'count', 'genetic_code': 27, 'max_steps': 10}, 200,
         'a peptide without instances (note no_stop_codon) of at most k bases is answered exact '
         '0 without a search, also after the request\'s budget stopped: stop null, work 0, '
         'determinism full, unlike the patterns a sticky stop leaves unknown (SPEC §7.6)',
         entries(field('stop', {'phase': 'discovery', 'reason': 'max_steps'}),
                 expect_all(exact(0), field('stop', None), has_note('no_stop_codon'),
                            work_zero))),
    post('genetic_code_unknown', 'masked',
         {'patterns': [p(NDM_PEP, 'protein')], 'genetic_code': 7}, 400,
         '400 genetic_code_unknown: genetic_code is not an NCBI translation table id (7 was '
         'merged into 4; capabilities genetic_codes lists the ids)',
         refused('genetic_code_unknown')),

    # ---------------------------------------------------------------- without a mask (#16)
    post('unmasked_count', 'unmasked_unchecked',
         {'patterns': [p(NDM_F, ident='NDM-F'), p(START16, ident='start'),
                       p(ABSENT_20, ident='absent')], 'mode': 'count'}, 200,
         'a graph without its dummy-edge mask (owner decision #16; it answered 400 mask_required '
         'before), the check of decision #24 off (--pattern-max-checked-entries 0: what a '
         'pattern above the limit answers): index counting upper_bound and dummy_fraction; a '
         'count the search could not resolve is bounds [lower, upper], upper the candidate '
         'entries (the source dummies among them), lower what was spelled whole (offset 0), '
         'each with its estimate round(upper x f) and the note estimate_sampled_dummy_fraction. '
         'NDM-F: [2, 24], estimate 24 (exact 24 with the mask); an island start: [2, 32], '
         'estimate 32, while the masked graph counts 17 -- the estimate is not a bound; an '
         'absent primer: an empty block, exact 0 (at the default limit the first two are '
         'checked: unmasked_checked)',
         entries(estimated(2, 24, 24), estimated(2, START16_UPPER, START16_UPPER),
                 expect_all(exact(0), field('notes', [])))),
    post('unmasked_checked', 'unmasked',
         {'patterns': [p(NDM_F, ident='NDM-F'), p(START16, ident='start'),
                       p(START14, ident='start14'), p(ABSENT_20, ident='absent')],
          'mode': 'count'}, 200,
         'owner decision #24 on the default server (max_checked_entries 50): a pattern with at '
         'most 50 unchecked candidates has each tested at query time (k - 1 = 30 steps each, in '
         'work.steps and mask_scans): NDM-F, [2, 24] unchecked (unmasked_count), its 22 '
         'candidates tested: exact 24, the masked count; the island start, [2, 32]: its 30 '
         'tested, the 15 source dummies dropped: exact 17, no estimate, no note; the first 14 '
         'bases of that start, [4, 72], 68 unchecked: above the limit, bounds with its estimate '
         'as with the check off; the absent primer exact 0',
         entries(expect_all(exact(24), field('notes', [])),
                 expect_all(exact(START16_EXACT), field('notes', [])),
                 estimated(START14_LOWER, START14_UPPER, START14_UPPER),
                 expect_all(exact(0), field('notes', [])))),
    post('unmasked_checked_dummies', 'unmasked',
         {'patterns': [p(START16, ident='start')], 'mode': 'count', 'scope': 'suffix',
          'strands': 'forward'}, 200,
         'owner decision #24: the island start in suffix scope on the forward strand has one '
         'candidate, the source dummy that holds it after 15 $ ([0, 1] with the check off); '
         'tested, it is a dummy: exact 0, an absence claim (absence_scope suffix_only) where the '
         'check off states only bounds',
         entries(expect_all(exact(0), field('notes', [])))),
    post('unmasked_labels_all', 'unmasked_unchecked',
         {'patterns': [p(NDM_F, ident='NDM-F'), p(ABSENT_20, ident='absent')],
          'output': {'labels': 'all'}}, 200,
         'all_or_count with labels "all" without the mask, the check of decision #24 off: the '
         'upper bound (24) within max_contexts, so every candidate is enumerated and the source '
         'dummies dropped: the list is exact (the masked graph\'s 24 contexts, 9 labels, 42 '
         'placed occurrences), and the counts, made from it, exact; the absent primer complete',
         entries(expect_all(exact(24), complete, labelled, counted('labels', 'exact', 9),
                            counted('occurrences', 'exact', 42)),
                 expect_all(exact(0), complete, counted('labels', 'exact', 0)))),
    post('unmasked_threshold_upper_bound', 'unmasked_unchecked',
         {'patterns': [p(START16, ident='start')], 'max_contexts': START16_EXACT,
          'output': {'labels': 'none'}}, 200,
         'all_or_count admits on the upper bound without the mask (conservative), the check of '
         'decision #24 off: the island start\'s 17 contexts would fit max_contexts 17 (the '
         'masked graph returns them), its upper bound 32 does not: withheld '
         'count_above_threshold with the count bounds, and the note threshold_upper_bound says '
         'that the lower bound was within the threshold (at the default limit its 30 unchecked '
         'candidates are tested and the 17 contexts released)',
         entries(expect_all(withheld('count_above_threshold'),
                            estimated(2, START16_UPPER, START16_UPPER),
                            has_note('threshold_upper_bound')))),
    post('unmasked_stop_at_threshold', 'unmasked',
         {'patterns': [p(NDM_F, ident='NDM-F')], 'max_contexts': 5, 'stop_at_threshold': True,
          'output': {'labels': 'none'}}, 200,
         'stop_at_threshold without the mask compares the running upper bound (also with the '
         'check of decision #24 on: the stop comes before it): discovery stops earlier than on '
         'the masked graph (at_least 0 here, 6 there), withheld threshold_crossed, note '
         'threshold_upper_bound (the lower bound had not crossed)',
         entries(expect_all(withheld('threshold_crossed'), relation('at_least'),
                            field('stop', {'phase': 'discovery', 'reason': 'max_contexts'}),
                            has_note('threshold_upper_bound')))),
    post('unmasked_partial', 'unmasked_unchecked',
         {'patterns': [p(NDM_F, ident='NDM-F')], 'mode': 'partial', 'max_contexts': 5,
          'output': {'labels': 'none'}}, 200,
         'partial without the mask, the check of decision #24 off: the first 5 contexts (the '
         'masked graph\'s first 5, no source dummy among them), cut max_contexts; the release '
         'raises each lower bound to what it released at that strand and offset: bounds [6, 24]',
         entries(expect_all(cut('max_contexts'), field('returned', 5), estimated(6, 24, 24)))),
    post('unmasked_paths', 'unmasked',
         {'patterns': [p(NDM_40, ident='NDM-40'), p(ABSENT_40, ident='absent')],
          'long_search': 'paths', 'output': {'labels': 'none'}}, 200,
         'patterns longer than k without the mask: an anchor window without a leading N is '
         'spelled whole, so the anchors are exact (no source dummy can hold one), and the paths '
         'are the masked graph\'s: exact, released, complete',
         entries(expect_all(exact(2), paths_count('exact', 2, 'completed'), complete),
                 expect_all(paths_count('exact', 0, 'no_anchors'), complete))),

    # ---------------------------------------------------------------- predicates
    post('predicate_filter', 'masked',
         {'patterns': [p(GCG12, ident='GCG12')], 'max_contexts': 100, 'predicate': SELECTIVE,
          'output': {'labels': 'predicate_only', 'occurrences': False}}, 200,
         'a selective predicate on a pattern above max_contexts (SPEC §19.13 a): the GCG repeat\'s '
         '1,828 contexts (withheld count_above_threshold without the predicate) tested in 1,777 '
         'rows, 60 selected and returned, complete; the predicate block (normal form, names, '
         'known, unknown_labels, vacuous, scope, strands either), selection {completed, kmer, '
         'rows}, absence_filter predicate; each result with selection_labels ["562"] and its own '
         'row\'s predicate labels (predicate_only, occurrences false: placement not_requested)',
         expect_all(predicate_block(normal_form=SELECTIVE, unknown=[], vacuous=False,
                                    strands='either'),
                    entries(expect_all(exact(1828), complete, field('returned', SELECTIVE_COUNT),
                                       selection('completed', ('exact', 1828),
                                                 ('exact', SELECTIVE_COUNT)),
                                       selection_labels_are(['562']),
                                       selection_strands_are(['both'], ['context']),
                                       lambda e: check(e['work']['predicate_rows'] == 1777,
                                                       e['work']))))),
    post('predicate_either', 'masked',
         {'patterns': [p(NDM_F, ident='NDM-F')], 'mode': 'count', 'predicate': {'none': ['546']},
          'output': {'labels': 'predicate_only'}}, 200,
         'P11, the strands of a BASIC graph: none(546) on the blaNDM-1 primer, predicate_strands '
         'either (the default): 546\'s records hold the primer on + only, so each of the 24 '
         'contexts has 546 on its k-mer or on its reverse complement: 0 selected; mode count '
         'reads the rows (24, 12 reverse-complement lookups), and states the projection it '
         'named and did not build (projection_not_read)',
         expect_all(predicate_block(vacuous=True, strands='either'),
                    entries(expect_all(exact(24), selection('completed', ('exact', 24),
                                                            ('exact', 0)),
                                       field('notes', ['projection_not_read']))))),
    post('predicate_context', 'masked',
         {'patterns': [p(NDM_F, ident='NDM-F')], 'mode': 'count', 'predicate': {'none': ['546']},
          'predicate_strands': 'context'}, 200,
         'the same with predicate_strands "context" (the deposited strand only): the 12 - '
         'contexts, whose own rows lack 546, are selected although 546\'s records carry the '
         'primer on the other strand',
         expect_all(predicate_block(strands='context'),
                    entries(selection('completed', ('exact', 24), ('exact', 12))))),
    post('predicate_at_least', 'masked',
         {'patterns': [p(NDM_F, ident='NDM-F')],
          'predicate': {'at_least': {'n': 8, 'labels': MINI_COLUMNS}},
          'predicate_strands': 'context'}, 200,
         'at_least 8 of the nine columns, "context": the 12 + contexts (9 columns on their rows; '
         'the - contexts have 7), returned and complete',
         entries(expect_all(complete, selection('completed', ('exact', 24), ('exact', 12)),
                            lambda e: check({r['strand'] for r in e['results']} == {'+'},
                                            e['results'])))),
    post('predicate_unknown_constant', 'masked',
         {'patterns': [p(NDM_F, ident='NDM-F')], 'predicate': {'any': [TYPO]}}, 200,
         'a typo (5622 is no column): reported in unknown_labels, never as an absence; the '
         'normal form is false (P18): no row read, selection constant, tested the 24 contexts, '
         'selected exact 0, complete; note predicate_constant',
         expect_all(predicate_block(normal_form=False, unknown=[TYPO]),
                    entries(expect_all(complete, field('returned', 0),
                                       selection('constant', ('exact', 24), ('exact', 0)),
                                       field('notes', ['predicate_constant']),
                                       lambda e: check(e['work']['predicate_rows'] == 0,
                                                       e['work']))))),
    post('predicate_unknown_folded', 'masked',
         {'patterns': [p(NDM_F, ident='NDM-F')],
          'predicate': {'and': [{'any': ['562']}, {'none': [TYPO]}]}}, 200,
         'an unknown name folded away (SPEC §19.4: none of an unknown is true): normal form '
         '{"any": ["562"]}, unknown_labels ["5622"]; the 24 contexts selected and returned',
         expect_all(predicate_block(normal_form={'any': ['562']}, unknown=[TYPO]),
                    entries(expect_all(complete, selection('completed', ('exact', 24),
                                                           ('exact', 24)))))),
    post('predicate_vacuous', 'masked',
         {'patterns': [p(GCG12, ident='GCG12')], 'mode': 'partial', 'max_contexts': 5,
          'predicate': {'none': ['287']}}, 200,
         'a vacuous predicate (true on a context carrying none of its labels: vacuous true), '
         'partial: 68 of the 1,828 contexts selected, the first 5 returned, cut max_contexts',
         expect_all(predicate_block(vacuous=True),
                    entries(expect_all(cut('max_contexts'), field('returned', 5),
                                       selection('completed', ('exact', 1828),
                                                 ('exact', 68)))))),
    post('predicate_budget', 'masked',
         {'patterns': [p(GCG12, ident='GCG12'), p(NDM_F, ident='NDM-F')], 'predicate': SELECTIVE,
          'predicate_strands': 'context', 'max_predicate_work': 1}, 200,
         'max_predicate_work 1, "context" (SPEC §19.13 d): the first row is read (the budget is '
         'checked before a read), the next is not: tested exact 1, selected bounds [0, 1827], '
         'stop {selection, max_predicate_work}, withheld predicate_budget; the stop is sticky: '
         'the next pattern\'s discovery runs, its selection not_started, predicate_budget',
         entries(expect_all(withheld('predicate_budget'), exact(1828),
                            selection('stopped', ('exact', 1), (0, 1827)),
                            field('stop', {'phase': 'selection',
                                           'reason': 'max_predicate_work'})),
                 expect_all(withheld('predicate_budget'), exact(24),
                            selection('not_started', ('unknown', None), ('unknown', None)),
                            field('stop', {'phase': 'selection',
                                           'reason': 'max_predicate_work'})))),
    post('predicate_budget_partial', 'masked',
         {'patterns': [p(GCG12, ident='GCG12')], 'mode': 'partial', 'predicate': SELECTIVE,
          'predicate_strands': 'context', 'max_predicate_work': 1}, 200,
         'the same in partial: nothing selected among the one context tested, returned 0, cut '
         'max_predicate_work',
         entries(expect_all(cut('max_predicate_work'), field('returned', 0),
                            selection('stopped', ('exact', 1), (0, 1827))))),
    post('predicate_above_threshold', 'masked',
         {'patterns': [p(GCG12, ident='GCG12')], 'predicate': SELECTIVE,
          'max_predicate_contexts': 100}, 200,
         'the compute admission: 1,828 raw contexts above max_predicate_contexts 100, '
         'all_or_count reads nothing: withheld predicate_above_threshold, selection '
         'not_admitted, tested and selected unknown, work.predicate_rows 0',
         entries(expect_all(withheld('predicate_above_threshold'),
                            selection('not_admitted', ('unknown', None), ('unknown', None)),
                            lambda e: check(e['work']['predicate_rows'] == 0, e['work'])))),
    post('predicate_partial_admission', 'masked',
         {'patterns': [p(GCG12, ident='GCG12')], 'mode': 'partial', 'predicate': SELECTIVE,
          'max_predicate_contexts': 100}, 200,
         'the same admission in partial: the first 100 raw contexts in answer order are tested '
         '(the count stays exact 1,828), selected bounds [S, S + 1,728], cut '
         'max_predicate_contexts',
         entries(expect_all(cut('max_predicate_contexts'), exact(1828),
                            lambda e: check(e['selection']['pass'] == 'stopped'
                                            and e['counts']['tested']['value'] == 100
                                            and e['counts']['selected']['upper']
                                            == e['counts']['selected']['lower'] + 1728,
                                            e['counts'])))),
    post('predicate_selected_above', 'masked',
         {'patterns': [p(GCG12, ident='GCG12')], 'max_contexts': 10, 'predicate': SELECTIVE},
         200,
         'more selected (60) than max_contexts (10): all_or_count withholds '
         'selected_above_threshold, the counts kept (selected exact 60)',
         entries(expect_all(withheld('selected_above_threshold'),
                            selection('completed', ('exact', 1828),
                                      ('exact', SELECTIVE_COUNT))))),
    post('predicate_stop_at_threshold', 'masked',
         {'patterns': [p(GCG12, ident='GCG12')], 'mode': 'partial', 'max_contexts': 10,
          'stop_at_threshold': True, 'predicate': SELECTIVE}, 200,
         'stop_at_threshold on the selected count, partial: the pass ends once 11 are '
         'selected (an early exit: 269 of the 1,777 rows read), stop {selection, max_contexts}, '
         'the first 10 returned, cut max_contexts, selected bounds',
         entries(expect_all(cut('max_contexts'), field('returned', 10),
                            field('stop', {'phase': 'selection', 'reason': 'max_contexts'}),
                            lambda e: check(e['selection']['pass'] == 'stopped'
                                            and e['counts']['selected']['relation'] == 'bounds',
                                            e['counts'])))),
    post('predicate_stop_at_raw', 'masked',
         {'patterns': [p(GCG12, ident='GCG12')], 'stop_at_threshold': True,
          'max_predicate_contexts': 100, 'predicate': SELECTIVE}, 200,
         'stop_at_threshold on the raw count with a predicate: discovery stops once its running '
         'count passes max_predicate_contexts (stop {discovery, max_predicate_contexts}), '
         'withheld threshold_crossed, the selection not_started',
         entries(expect_all(withheld('threshold_crossed'), relation('at_least'),
                            field('stop', {'phase': 'discovery',
                                           'reason': 'max_predicate_contexts'}),
                            selection('not_started', ('unknown', None), ('unknown', None))))),
    post('predicate_only_record', 'masked',
         {'patterns': [p(NDM_F, ident='NDM-F')], 'predicate': {'any': ['562']},
          'output': {'labels': 'predicate_only'}}, 200,
         'output.labels "predicate_only" with occurrences (the default): each selected context '
         'with the predicate\'s labels on its own row, placed in the records (562\'s 13 placed '
         'occurrences over the 24 contexts, as labels_all\'s by_label has them), by_label of '
         'those labels only',
         entries(expect_all(complete, selection('completed', ('exact', 24), ('exact', 24)),
                            counted('labels', 'exact', 1), counted('occurrences', 'exact', 13),
                            field('placement', 'record'), selection_labels_are(['562']),
                            selection_strands_are(['both'])))),
    post('predicate_selection_strands', 'masked',
         {'patterns': [p(NDM_F, ident='NDM-F')], 'predicate': {'any': ['546']},
          'output': {'labels': 'predicate_only', 'occurrences': False}}, 200,
         'per selected result and label the orientation whose row carries it (the owner\'s answer '
         'to P11, selection_strands): 546 holds the blaNDM-1 primer on + only, so any(546) '
         '"either" selects all 24 contexts, the 12 + ones by their own row ("context", their '
         'labels [546]) and the 12 - ones by their reverse complement\'s ("reverse_complement", '
         'their own labels empty)',
         entries(expect_all(complete, selection('completed', ('exact', 24), ('exact', 24)),
                            selection_labels_are(['546']),
                            selection_strands_are(['context'], ['reverse_complement']),
                            lambda e: check(all(
                                r['selection_strands'] == (['context'] if r['strand'] == '+'
                                                           else ['reverse_complement'])
                                and len(r['labels']) == (r['strand'] == '+')
                                for r in e['results']), e['results'])))),
    post('predicate_all', 'masked',
         {'patterns': [p(NDM_F, ident='NDM-F')], 'predicate': {'none': ['546']},
          'predicate_strands': 'context', 'output': {'labels': 'all'}}, 200,
         'a predicate with output.labels "all" (P8): the 12 - contexts selected by none(546) '
         '"context", each with every label of its row (7), placed; their selection_labels empty '
         '(none of the predicate\'s labels is on them)',
         entries(expect_all(complete, selection('completed', ('exact', 24), ('exact', 12)),
                            counted('labels', 'exact', 7), selection_labels_are([]),
                            selection_strands_are([]),
                            lambda e: check(all(r['labels_total'] == 7 for r in e['results']),
                                            e['results'])))),
    post('predicate_long_anchors', 'masked',
         {'patterns': [p(NDM_40, ident='NDM-40')], 'predicate': {'any': ['562']}}, 200,
         'a pattern longer than k with a predicate under long_search "anchors" (the default): '
         'its anchors\' answer as without it (withheld paths_later_increment), its selection '
         'not_started (a predicate selects supported paths, a later increment)',
         entries(expect_all(exact(2), withheld('paths_later_increment'),
                            lambda e: check(e['selection']['pass'] == 'not_started'
                                            and e['counts']['selected']['relation']
                                            == 'unknown', e)))),
    post('predicate_unmasked', 'unmasked_unchecked',
         {'patterns': [p(NDM_F, ident='NDM-F')], 'predicate': {'any': ['562']}}, 200,
         'a predicate on the graph without its mask, the check of decision #24 off: the raw '
         'count bounds [2, 24] is admitted on its upper bound, the release enumerates every '
         'candidate and drops the source dummies, so the raw count is exact 24; the 24 '
         'contexts selected and returned, complete',
         entries(expect_all(exact(24), complete,
                            selection('completed', ('exact', 24), ('exact', 24))))),
    post('predicate_unbudgeted', 'primary',
         {'patterns': [p(NDM_F)], 'mode': 'count', 'predicate': {'any': ['562.fa']}}, 400,
         '400 annotation_unbudgeted: a predicate reads the annotation in every mode, and the '
         'PRIMARY index\'s column annotation has no budget-aware decode',
         refused('annotation_unbudgeted')),
    post('predicate_unbudgeted_allowed', 'primary',
         {'patterns': [p(NDM_F, ident='NDM-F')], 'predicate': {'none': ['573.fa']},
          'predicate_strands': 'context', 'allow_unbudgeted_annotation': True}, 200,
         'the same with allow_unbudgeted_annotation: one row serves a k-mer and its reverse '
         'complement on a PRIMARY graph, so the predicate block says strands either whatever '
         'was asked (limits.predicate_strands: context, as requested); selection.access '
         'columns (single cells for at most 16 labels); note annotation_unbudgeted',
         expect_all(predicate_block(strands='either'),
                    entries(expect_all(lambda e: check(e['selection']['access'] == 'columns',
                                                       e['selection']),
                                       has_note('annotation_unbudgeted'))),
                    lambda a: check(a['limits']['predicate_strands'] == 'context',
                                    a['limits']))),
    post('predicate_invalid', 'masked',
         {'patterns': [p(NDM_F)], 'predicate': {'any': [562]}}, 400,
         '400 invalid_request: a name is a string (P7; the message names the path and the '
         'fix, "write a taxid as \\"562\\"")',
         refused('invalid_request')),
    post('predicate_paths_refused', 'masked',
         {'patterns': [p(NDM_40)], 'predicate': {'any': ['562']}, 'long_search': 'paths'}, 400,
         '400 invalid_request: a predicate with long_search "paths" (P24: a predicate selects '
         'among supported paths, long_search "supported_paths", a later increment)',
         refused('invalid_request')),
    post('predicate_too_large', 'masked_small_predicate_cap',
         {'patterns': [p(NDM_F)], 'predicate': {'any': MINI_COLUMNS[:SMALL_PREDICATE_CAP + 1]}},
         400,
         '400 predicate_too_large on a server whose predicates may list 4 names '
         '(--pattern-max-predicate-labels 4, capabilities caps.max_predicate_labels): 5 names; '
         'the message names the count and the cap',
         refused('predicate_too_large')),

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
    post('predicate_only_without_predicate', 'masked',
         {'patterns': [p(NDM_F)], 'output': {'labels': 'predicate_only'}}, 400,
         '400 invalid_request: output.labels "predicate_only" returns the labels a predicate '
         'names, so it needs one (served with a predicate since increment 5b; 400 '
         'later_increment before, fixture later_increment_labels)',
         refused('invalid_request')),
    post('later_increment_graphs', 'masked',
         {'patterns': [p(NDM_F)], 'mode': 'count', 'graphs': [INDEX_NAME]}, 400,
         '400 later_increment: the request field graphs (multi-graph selection) is refused by '
         'name, whatever its value',
         refused('later_increment')),
    post('resident_only', 'masked',
         {'patterns': [p(NDM_F)], 'mode': 'count', 'in_ram': False}, 400,
         '400 resident_only: in_ram, whatever its value (the route never loads an index)',
         refused('resident_only')),
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
                + '  (the mini index as built: no .edgemask; counting upper_bound; at most 50 '
                  'unchecked candidates tested, the default)',
    'unmasked_unchecked': 'server_query -i {mini}/graph_k31.dbg -a {mini}/' + MINI_ANNO
                          + ' --pattern-max-checked-entries 0 --index-name ' + INDEX_NAME
                          + ' --index-release ' + INDEX_RELEASE
                          + '  (the same, no unchecked candidate tested: every count with '
                            'unchecked candidates bounds)',
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
    'masked_small_predicate_cap': 'server_query -i {work}/masked/graph_k31.dbg -a {work}/masked/'
                                  + MINI_ANNO + ' --pattern-max-predicate-labels '
                                  + str(SMALL_PREDICATE_CAP) + ' --index-name ' + INDEX_NAME
                                  + ' --index-release ' + INDEX_RELEASE
                                  + '  (the masked copy; a predicate may list '
                                  + str(SMALL_PREDICATE_CAP) + ' names)',
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
    if name == 'unmasked_unchecked':
        return ['-i', os.path.join(mini, MINI_GRAPH), '-a', os.path.join(mini, MINI_ANNO),
                '--pattern-max-checked-entries', '0'] + ident
    if name == 'built_at_load':
        return ['-i', os.path.join(mini, MINI_GRAPH), '-a', os.path.join(mini, MINI_ANNO),
                '--pattern-build-mask'] + ident
    if name == 'multi':
        # one name and one manifest describe one index: a multi-graph server takes neither
        return [os.path.join(work, 'graphs.csv'), '--index-release', INDEX_RELEASE]
    if name == 'masked_no_map':
        return ['-i', os.path.join(masked, MINI_GRAPH), '-a', os.path.join(masked, MINI_ANNO),
                '--no-coord-mapping'] + ident
    if name == 'masked_small_predicate_cap':
        return ['-i', os.path.join(masked, MINI_GRAPH), '-a', os.path.join(masked, MINI_ANNO),
                '--pattern-max-predicate-labels', str(SMALL_PREDICATE_CAP)] + ident
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
                # the same in every run
                clocked = True
                e['counts'] = blank_count(e['counts'])
                e['work'] = {k: BLANK for k in e['work']}
                if e.get('cut') == {'reason': 'time'} \
                        and e.get('stop') != {'phase': 'discovery', 'reason': 'time'}:
                    # the clock cut its release: how many contexts got out (a time stop in
                    # discovery releases nothing, returned 0, compared)
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
        '`mask_required` is retired (owner decision #16 of 2026-10-08): a graph without its',
        'dummy-edge mask is served with upper bounds and estimates (the `unmasked_*` fixtures),',
        'and no build since answers it; test_pattern_fixtures.py names it in RETIRED. A pattern',
        'with few unchecked candidates there (at most `caps.max_checked_entries`, 50 by default)',
        'has each tested and its counts exact (owner decision #24: `unmasked_checked`,',
        '`unmasked_checked_dummies`); the server `unmasked_unchecked` has the check off.',
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
