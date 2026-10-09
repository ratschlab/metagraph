"""MGT v1 reader and canonical writer (DESIGN-traverse-graphlet.md §2).

parse(text) reads one document in one pass (records in document order), then resolves
each arm's derived fields in segment-id order once its R records are known (G.end '*'
needs the runs anchored on the segment; a child's G.entry delta needs its parent's end
set). It validates the A counts and the Z line count, so a truncated body raises
GraphletFormatError instead of producing a shallower graphlet (§14: a cut-off parse is a
local failure, never a result).

dump(g) writes the canonical form (§2.3): '*' wherever the field's rule reproduces
the value, every SETEXPR in its shortest form (tie explicit), RANGES collapsed, the
B/T maps in canonical order. dump(parse(x)) == x for every canonical document;
is_canonical(x) tests exactly that.
"""

import hashlib
import json
import os
import re
import tempfile
from array import array

from . import budget as _B
from . import coords, derive
from ._codec import (
    CodecError, GraphletFormatError, LabelSetInterner, QUAL, REASON, REASON_ORDER,
    RESOURCE_CODES, decode_kvalue, decode_pairs, decode_ranges, decode_setexpr,
    encode_kvalue, encode_pairs, encode_ranges, encode_setexpr, end_token, fmt_bool,
    fmt_int, fmt_num, front_decode_parts, front_encode, parse_bool, parse_end_token,
    parse_int, parse_num, pct_escape, pct_unescape, tok as _tok,
)
from .budget import (ARRAY, LIST, LIST_ITEM, STR, TEXT_LINE, W_DUMP_LINE, W_LINE,
                     LocalBudgetExceeded)
from .model import (
    Arm, BranchEvent, CapTrigger, Counts, Dropped, Event, Frontier, Graphlet, GrowthBin,
    Label, LabelsPerNode, Leaf, LeafContinuation, Limitation, Outcome, PresenceRun,
    Refusal, ResourceStop, Run, Segment, SeedInfo,
)

# local limits: the account of a parse, fitted to the traced peak of the corpus (CPython 3.11,
# 140 bodies; actual/fit at most 1.47, 1.25 on peaks of 100 KB or more) and charged at 1.5
# times the fit: per body byte and per line (the line list, the tokens), per record of the
# types that hold the most (a segment's raw tokens until its arm resolves, a leaf's
# extras, the K and X records' typed values), and per label set the parse interns for the
# first time (the array and its key)
_MODEL_PER_BYTE = 3.1
_MODEL_PER_LINE = 540
_MODEL_PER_RECORD = (('\nG ', 950), ('\nT ', 640), ('\nR ', 65), ('\nV ', 71),
                     ('\nK ', 3950), ('\nX ', 2390))
_MODEL_PER_SET = 124
_MODEL_PER_SET_ID = 13.3
_LINE_BLOCK = 64

__all__ = ['parse', 'dump', 'is_canonical', 'load', 'save', 'from_response',
           'standalone_text', 'seed_envelope', 'RESERVED_ENVELOPE_NAMES',
           'GraphletFormatError', 'FORMAT_VERSION', 'utf8_bytes']

FORMAT_VERSION = 1

_MODE = {'c': 'constrain', 'a': 'annotate'}
_SUPPORT = {'k': 'kmer', 't': 'trace'}
_RECONVERGE = {'m': 'merge', 'k': 'keep'}
_STATUS = {'c': 'complete', 't': 'truncated', 'p': 'pruned'}
_SCOPE = {'p': 'per_path', 'u': 'united_history'}
_KIND = {'c': 'column', 'h': 'header'}
_O_WALKS = {'c': 'complete', 'p': 'partial', 'f': 'failed'}
_O_BRANCH = {'c': 'complete', 'x': 'cut'}
# label_evidence 'qualified' is 'q' (the O grammar of §2.2 does not list it)
_O_LABEL = {'c': 'complete', 'l': 'lower_bound', 'q': 'qualified'}
_O_DELIVERY = {'i': 'inline', 's': 'spooled', 'p': 'paged'}
_REGIMES = ('basic', 'primary', 'canonical')
_TOKEN_RE = re.compile(r'[!-~]+\Z')          # printable ASCII without space
_NS_RE = re.compile(r'[A-Za-z0-9._-]+\Z')
_NAME_RE = re.compile(r'[A-Za-z0-9_.\[\]-]+\Z')   # counter / extra names


# which fields a record wrote as '*' (check_rules): one shared frozenset per combination
# -- a new frozenset per run or segment would cost 216 B each (frozenset() is no
# singleton), 3.4 MB of a 13.9 MB model on 15,922 runs
_NO_STARS = frozenset()
_STRUCT_STAR = frozenset(['structural_successors'])
_STAR_SETS = {}


def _star_set(names):
    key = frozenset(names)
    return _STAR_SETS.setdefault(key, key)


def _inv(d):
    return {v: k for k, v in d.items()}


def _code(table, tok, what):
    try:
        return table[tok]
    except KeyError:
        raise CodecError('%s %s' % (what, _tok(tok))) from None


def _opt(tok, conv):
    return None if tok == '*' else conv(tok)


def _star(v, conv):
    return '*' if v is None else conv(v)


def _reason_code(tok):
    if len(tok) != 1 or tok not in REASON:
        raise CodecError('end-reason code %s' % _tok(tok))
    return tok


def _csv_ints(tok):
    if tok == '.':
        return []
    return [parse_int(x) for x in tok.split(',')]


def _fixed(tok, what):
    if not _TOKEN_RE.match(tok):
        raise CodecError('%s %s is not printable ASCII without spaces' % (what, _tok(tok)))
    return tok


# ======================================================================= reading

class _ArmRaw:
    """An arm's records as read, before the derived fields are resolved."""
    __slots__ = ('arm', 'seg_raw', 'runs', 'run_lines')

    def __init__(self, arm):
        self.arm = arm
        self.seg_raw = []      # per segment: dict of raw tokens + line numbers
        self.runs = []
        self.run_lines = []    # the line of each R record (the checks after the read)


class _Reader:
    def __init__(self, lines, b=None):
        self.lines = lines
        self.i = 0
        self.interner = LabelSetInterner()
        # local limits: lines are paid for in blocks admitted before they are read
        self.b = b
        self.paid = 0
        self.per_line_work = 0
        self.per_line_bytes = 0

    def _pay(self):
        n = min(_LINE_BLOCK, len(self.lines) - self.i)
        self.b.done['lines'] = self.i
        # the block's share of the work and of the modelled model (spread evenly)
        self.b.charge((W_LINE + self.per_line_work) * n, int(self.per_line_bytes * n))
        self.paid = self.i + n

    # ------------------------------------------------------------ line access

    def peek_tag(self):
        if self.i >= len(self.lines):
            return None
        line = self.lines[self.i]
        return line[:1] if len(line) >= 2 and line[1] == ' ' else line

    def take(self, tag, n=None, maxsplit=None):
        """The next line's fields; |n| exact count (with |maxsplit| the free-text tail
        is the last field)."""
        if self.i >= len(self.lines):
            raise GraphletFormatError(0, 'truncated document: expected %s, found the end'
                                      % tag)
        if self.b is not None and self.i >= self.paid:
            self._pay()
        line = self.lines[self.i]
        self.i += 1
        if maxsplit is not None:
            parts = line.split(' ', maxsplit)
        else:
            parts = line.split(' ')
        if parts[0] != tag:
            raise GraphletFormatError(self.i, 'expected a %s record, found %s'
                                      % (tag, _tok(line)))
        if n is not None and (len(parts) != n if not isinstance(n, tuple)
                              else len(parts) not in n):
            raise GraphletFormatError(self.i, '%s record with %d fields' % (tag, len(parts)))
        if maxsplit is None and '' in parts:
            raise GraphletFormatError(self.i, 'empty field (two spaces) in a %s record' % tag)
        return parts

    def err(self, msg, line_no=None):
        return GraphletFormatError(self.i if line_no is None else line_no, msg)

    # ------------------------------------------------------------ the document

    def run(self):
        try:
            return self._run()
        except CodecError as e:
            raise GraphletFormatError(self.i, str(e)) from None

    def _run(self):
        h = self.take('H', 16)
        if h[1] != 'mgt':
            raise self.err('not an MGT document')
        if h[2] != str(FORMAT_VERSION):
            raise self.err('unsupported MGT version %s (this reader: %d)'
                           % (_tok(h[2]), FORMAT_VERSION))
        if h[4] not in _REGIMES:
            raise self.err('regime %s' % _tok(h[4]))
        if h[12] != 'walk':
            raise self.err('orientation %s (MGT v1 stores walking order: "walk")' % _tok(h[12]))
        if h[13] != '*' and not _NS_RE.match(h[13]):
            raise self.err('index_ns %s' % _tok(h[13]))
        g = Graphlet(
            format=FORMAT_VERSION, k=parse_int(h[3]), regime=h[4],
            alphabet=_fixed(h[5], 'alphabet'), mode=_code(_MODE, h[6], 'mode'),
            support=_code(_SUPPORT, h[7], 'support'),
            reconverge=_code(_RECONVERGE, h[8], 'reconverge'), cap=parse_int(h[9]),
            continuation_bp=parse_int(h[10]), seed_index=parse_int(h[11]),
            index_ns=_opt(h[13], str), index_fp=_opt(h[14], str),
            index_meta_fp=_opt(h[15], str), seed=None, labels=[], dropped=[],
            outcome=None, resource_stop=None, limitations=[], arms={})
        if self.peek_tag() == 'J':
            if self.b is not None:
                if self.i >= self.paid:
                    self._pay()
                # the envelope as JSON objects: a few times its text
                self.b.charge(len(self.lines[self.i]) >> 4, 8 * len(self.lines[self.i]))
            line = self.lines[self.i]
            self.i += 1
            try:
                j = json.loads(line[2:])
            except ValueError as e:
                raise self.err('J line is not JSON: %s' % e) from None
            if not isinstance(j, dict):
                raise self.err('J line is not a JSON object')
            _attach_j(g, j)
            g.has_j = True

        s = self.take('S', 6)
        g.seed = SeedInfo('' if s[1] == '*' else _fixed(s[1], 'seed id'), parse_int(s[2]),
                          parse_int(s[3]), parse_int(s[4]), _fixed(s[5], 'seed sequence'))
        if len(g.seed.sequence) != g.seed.length_bp:
            raise self.err('S: the sequence has %d bases, length_bp says %d'
                           % (len(g.seed.sequence), g.seed.length_bp))

        while self.peek_tag() == 'X':
            x = self.take('X', 4, maxsplit=3)
            g.dropped.append(Dropped(_fixed(x[1], 'reason'), decode_pairs(x[2]),
                                     pct_unescape(x[3])))

        prev = ''
        while self.peek_tag() == 'L':
            f = self.take('L', 6, maxsplit=5)
            kind = _code(_KIND, f[1], 'label kind')
            column = parse_int(f[2])
            if kind == 'header':
                seq_id = parse_int(f[3])
            elif f[3] != '*':
                raise self.err('a column label has no seq_id (write *)')
            else:
                seq_id = None
            if f[4] == '' or not re.match(r'(?:0|[1-9][0-9]*)\Z', f[4]):
                raise self.err('L prefix_len %s' % _tok(f[4]))
            name = front_decode_parts(prev, f[4], f[5])
            g.labels.append(Label(len(g.labels), kind, column, seq_id, name))
            prev = name
        if g.mode == 'constrain' and g.seed.num_seed_labels > len(g.labels):
            raise self.err('S: %d seed labels but %d L records'
                           % (g.seed.num_seed_labels, len(g.labels)))

        o = self.take('O', 5)
        g.outcome = Outcome(_code(_O_WALKS, o[1], 'outcome walks'),
                            _code(_O_BRANCH, o[2], 'outcome branch_diagnostics'),
                            _code(_O_LABEL, o[3], 'outcome label_evidence'),
                            _code(_O_DELIVERY, o[4], 'outcome delivery'))
        if self.peek_tag() == 'Q':
            q = self.take('Q', 10, maxsplit=9)
            for k in range(1, 9):
                if not _TOKEN_RE.match(q[k]):
                    raise self.err('Q field %d %s' % (k, _tok(q[k])))
            g.resource_stop = ResourceStop(
                q[1], q[2], q[3], *(_opt(t, decode_kvalue) for t in q[4:8]),
                [] if q[8] == '.' else q[8].split(','), pct_unescape(q[9]))
        while self.peek_tag() == 'K':
            g.limitations.append(self._limitation())

        sides_seen = []
        while self.peek_tag() == 'A':
            arm = self._arm(g)
            if arm.side in sides_seen or (sides_seen and arm.side == 'left'):
                raise self.err('arms out of order (left before right, once each)')
            sides_seen.append(arm.side)
            g.arms[arm.side] = arm

        z = self.take('Z', 2)
        if self.i != len(self.lines):
            raise GraphletFormatError(self.i + 1, 'records after Z')
        if parse_int(z[1]) != len(self.lines):
            raise GraphletFormatError(self.i, 'Z says %s lines, the document has %d'
                                      % (z[1], len(self.lines)))
        return g

    def _limitation(self):
        f = self.take('K', 9, maxsplit=8)
        arm = {'l': 'left', 'r': 'right', '*': None}.get(f[1], False)
        if arm is False:
            raise self.err('K arm %s' % _tok(f[1]))
        for k in (2, 3):
            _fixed(f[k], 'K field')
        extra = []
        if f[7] != '.':
            for item in f[7].split(','):
                name, eq, val = item.partition('=')
                if not eq or not _NAME_RE.match(name):
                    raise self.err('K extra item %s' % _tok(item))
                extra.append((name, decode_kvalue(val)))
        return Limitation(arm, f[2], f[3], decode_kvalue(f[4]), decode_kvalue(f[5]),
                          _opt(f[6], parse_int), extra, pct_unescape(f[8]))

    # ------------------------------------------------------------ one arm

    def _arm(self, g):
        a = self.take('A', 20)
        side = {'l': 'left', 'r': 'right'}.get(a[1])
        if side is None:
            raise self.err('A arm %s' % _tok(a[1]))
        counters = {}
        if a[10] != '.':
            for item in a[10].split(','):
                name, eq, val = item.partition('=')
                if not eq or not _NAME_RE.match(name) or name in counters:
                    raise self.err('A counter %s' % _tok(item))
                counters[name] = parse_int(val)
        cap_trigger = None
        if a[11] != '*':
            c = a[11].split(',')
            if len(c) != 7:
                raise self.err('A cap_trigger %s' % _tok(a[11]))
            cap_trigger = CapTrigger(_reason_code(c[0]), parse_int(c[1]), parse_int(c[2]),
                                     parse_int(c[3]), parse_int(c[4]), parse_bool(c[5]),
                                     parse_num(c[6]))
        arm = Arm(side=side, status=_code(_STATUS, a[2], 'status'),
                  complete_to_bp=parse_int(a[3]), scope=_code(_SCOPE, a[4], 'scope'),
                  frontier=Frontier(parse_int(a[5]), parse_int(a[6]), parse_bool(a[7])),
                  labels_per_node=LabelsPerNode(parse_int(a[8]), parse_int(a[9])),
                  cap_trigger=cap_trigger, counters=counters, growth=[], branch_events=[],
                  branch_events_total=parse_int(a[12]),
                  evidence_complete_to_bp=_opt(a[13], parse_int),
                  counts=Counts(*(parse_int(t) for t in a[14:20])))
        a_line = self.i
        while self.peek_tag() == 'B':
            arm.growth.append(self._growth())
        while self.peek_tag() == 'V':
            arm.branch_events.append(self._branch_event(len(g.labels)))
        raw = _ArmRaw(arm)
        while self.peek_tag() == 'G':
            raw.seg_raw.append(self._segment_records(g, len(raw.seg_raw)))
        while self.peek_tag() == 'R':
            raw.run_lines.append(self.i + 1)
            raw.runs.append(self._run_record(g, len(raw.runs), len(raw.seg_raw)))
        _resolve_arm(self, g, raw)
        c = arm.counts
        got = Counts(len(arm.segments), len(arm.runs), len(derive.leaves(arm)),
                     len(derive.splits(arm, g.mode)),
                     sum(1 for s in arm.segments if len(s.parents) > 1),
                     sum(s.length_bp for s in arm.segments))
        if got != c:
            raise GraphletFormatError(a_line, 'the %s arm does not match its A counts '
                                      '(truncated or corrupt body): A says %s, the body '
                                      'has %s' % (side, c, got))
        _check_arm_invariants(g, arm, a_line)
        return arm

    def _growth(self):
        b = self.take('B', 15)
        ends = {}
        if b[14] != '.':
            for item in b[14].split(','):
                code, colon, n = item.partition(':')
                if not colon or code in ends:
                    raise self.err('B label_ends item %s' % _tok(item))
                ends[_reason_code(code)] = parse_int(n)
        v = [parse_int(t) for t in b[1:5]]
        return GrowthBin(v[0], v[1], v[2], v[3], parse_bool(b[5]),
                         *(parse_int(t) for t in b[6:14]), ends)

    def _branch_event(self, n_labels):
        f = self.take('V', 8)
        refused = []
        if f[7] != '.':
            for item in f[7].split(';'):
                if len(item) < 3 or item[1] != ':':
                    raise self.err('V refusal %s' % _tok(item))
                cause, colon, ranges = item[2:].partition(':')
                if not colon or not cause:
                    raise self.err('V refusal %s' % _tok(item))
                refused.append(Refusal(item[0], cause, self._ids(decode_ranges(ranges, n_labels),
                                                                 n_labels)))
        chars = '' if f[3] == '.' else f[3]
        return BranchEvent(parse_int(f[1]), parse_int(f[2]), chars, _csv_ints(f[4]),
                           self._ids(decode_ranges(f[5], n_labels), n_labels),
                           self._ids(decode_ranges(f[6], n_labels), n_labels), refused)

    def _ids(self, ids, n_labels):
        if ids and ids[-1] >= n_labels:
            # a CodecError, so that the caller names the record's line (the resolution
            # pass runs after the arm has been read)
            raise CodecError('label id %d beyond the %d L records' % (ids[-1], n_labels))
        if self.b is None:
            return self.interner.intern(ids)
        # a decoded set is charged as it is interned: the decode (a delta copies its
        # whole base) and, for a set seen for the first time, its array
        misses = self.interner.misses
        self.b.charge(len(ids) >> _B.ID_SHIFT_PARSE)
        got = self.interner.intern(ids)
        if self.interner.misses != misses:
            self.b.charge(0, int(_MODEL_PER_SET + _MODEL_PER_SET_ID * len(ids)))
        return got

    def _segment_records(self, g, seg_id):
        line = self.i + 1
        f = self.take('G', 11)
        raw = {'line': line, 'g': f, 'P': [], 'E': [], 'T': None, 'C': None}
        while True:
            tag = self.peek_tag()
            if tag == 'P':
                if g.mode != 'annotate':
                    raise GraphletFormatError(self.i + 1, 'P records exist in annotate mode only')
                if raw['E'] or raw['T']:
                    raise GraphletFormatError(self.i + 1, 'P after E/T')
                raw['P'].append((self.i + 1, self.take('P', 5)))
            elif tag == 'E':
                if raw['T']:
                    raise GraphletFormatError(self.i + 1, 'E after T')
                raw['E'].append(self._event(len(g.labels)))
            elif tag == 'T':
                if raw['T']:
                    raise GraphletFormatError(self.i + 1, 'two T records in one segment')
                raw['T'] = (self.i + 1, self.take('T', 3))
            elif tag == 'C':
                if not raw['T'] or raw['C']:
                    raise GraphletFormatError(self.i + 1, 'C without its T record')
                raw['C'] = (self.i + 1, self.take('C', (5, 6)))
            else:
                return raw

    def _event(self, n_labels):
        line = self.lines[self.i]
        f = line.split(' ')
        if len(f) < 3:
            raise GraphletFormatError(self.i + 1, 'E record %s' % _tok(line))
        t = f[2]
        n = {'s': 6, 'b': 7, 'h': 7, 'v': 5, 't': 5, 'u': 5}.get(t)
        if n is None:
            raise GraphletFormatError(self.i + 1, 'event type %s' % _tok(t))
        f = self.take('E', n)
        at = parse_int(f[1])
        if t == 's':
            return Event(at, 'switch', from_label=self._label_id(f[3], n_labels),
                         to_label=self._label_id(f[4], n_labels), cost=parse_num(f[5]))
        if t == 'b':
            if len(f[3]) != 1:
                raise self.err('blocked char %s' % _tok(f[3]))
            return Event(at, 'blocked', char=f[3], reason=_reason_code(f[4]),
                         total=parse_int(f[5]),
                         labels=self._ids(decode_ranges(f[6], n_labels), n_labels))
        if t == 'h':
            if len(f[3]) != 1 or f[6] not in ('f', 's'):
                raise self.err('hairpin event %s' % _tok(line))
            return Event(at, 'hairpin', char=f[3], total=parse_int(f[4]),
                         labels=self._ids(decode_ranges(f[5], n_labels), n_labels),
                         followed=f[6] == 'f')
        if t == 'v':
            return Event(at, 'revisit', segment=parse_int(f[3]),
                         delta=None if f[4] == '=' else parse_int(f[4]))
        if t == 't':
            return Event(at, 'tip', char=f[3], length_bp=parse_int(f[4]))
        return Event(at, 'bubble', length_bp=parse_int(f[3]), alleles=_fixed(f[4], 'alleles'))

    def _label_id(self, tok, n_labels):
        i = parse_int(tok)
        if i >= n_labels:
            raise self.err('label id %d beyond the %d L records' % (i, n_labels))
        return i

    def _run_record(self, g, run_id, n_segments):
        f = self.take('R', (12, 13))
        code, qual, silent, merged = parse_end_token(f[5])
        if (len(f) == 13) != (code == 'B'):
            raise self.err('R needed_budget is present exactly for a loss_budget end')
        seg = parse_int(f[1])
        if seg >= n_segments:
            raise self.err('R anchor %d beyond the arm\'s %d segments' % (seg, n_segments))
        from_label = cost = None
        if f[7] != '*':
            fl, colon, c = f[7].partition(':')
            if not colon:
                raise self.err('R from_label:cost %s' % _tok(f[7]))
            from_label = self._label_id(fl, len(g.labels))
            cost = parse_num(c)
        prev = _opt(f[8], parse_int)
        if prev is not None and prev >= run_id:
            raise self.err('R prev_run %d is not an earlier run' % prev)
        if silent or merged:
            if f[9] != '*':
                raise self.err('R structural_successors is * for Lw / m runs')
            structural = None
        else:
            structural = parse_int(f[9])
        stars = _STRUCT_STAR if f[9] == '*' else _NO_STARS
        return Run(id=run_id, segment=seg, label=self._label_id(f[2], len(g.labels)),
                   from_bp=parse_int(f[3]), to_bp=parse_int(f[4]), end=f[5],
                   reason=None if merged else REASON[code], qualifier=qual, silent=silent,
                   merged=merged, route_bp=parse_int(f[6]), from_label=from_label, cost=cost,
                   prev_run=prev, structural_successors=structural,
                   branches=parse_int(f[10]), loss=parse_num(f[11]),
                   needed_budget=parse_num(f[12]) if len(f) == 13 else None, stars=stars)


def _resolve_arm(rd, g, raw):
    """Second pass over one arm: segment-id order (parents precede children), so every
    SETEXPR base and every '*' rule input is known when it is needed."""
    arm = raw.arm
    arm.runs = raw.runs
    n = len(raw.seg_raw)
    anchored = [[] for _ in range(n)]
    for r in arm.runs:
        if r.to_bp < r.from_bp:
            raise GraphletFormatError(raw.run_lines[r.id], '%s arm: run %d ends before it '
                                      'starts' % (arm.side, r.id))
        anchored[r.segment].append(r.id)
    arm.cache['runs_by_segment'] = anchored
    segs = arm.segments
    sequences = None
    for sid, sr in enumerate(raw.seg_raw):
        line = sr['line']
        f = sr['g']
        try:
            parents = () if f[1] == '*' else tuple(parse_int(t) for t in f[1].split(','))
            if any(p >= sid for p in parents) or len(set(parents)) != len(parents):
                raise CodecError('parents %s must be distinct earlier segments' % _tok(f[1]))
            if not parents and sid != 0:
                raise CodecError('only segment 0 is a root')
            if sid == 0 and parents:
                raise CodecError('segment 0 is the root')
            from_bp = parse_int(f[2])
            length = parse_int(f[3])
            for p in parents:
                if segs[p].end_bp != from_bp:
                    raise CodecError('segment starts at %d, its parent %d ends at %d'
                                     % (from_bp, p, segs[p].end_bp))
            if not parents and from_bp != 0:
                raise CodecError('the root starts at 0')
            partition = derive.rule_partition(parents) if f[7] == '*' else None
            if partition is None:
                parts = f[7].split('|')
                if len(parents) <= 1 or len(parts) != len(parents):
                    raise CodecError('partition %s for %d parent(s)' % (_tok(f[7]), len(parents)))
                partition = [rd._ids(decode_ranges(p, len(g.labels)), len(g.labels))
                             for p in parts]
            seg = Segment(id=sid, parents=parents, from_bp=from_bp, length_bp=length,
                          walk=None, first_base=None, entry=None, entry_total=0, end=None,
                          partition=partition, split=_opt(f[8], parse_bool))
            stars = set()
            for name, tok in (('parents', f[1]), ('entry', f[4]), ('entry_total', f[5]),
                              ('end', f[6]), ('partition', f[7])):
                if tok == '*':
                    stars.add(name)
            seg.stars = _star_set(stars)
            if f[4] == '*':
                rule = derive.rule_entry(g.mode, g.seed.num_seed_labels, seg)
                if rule is None:
                    raise CodecError('G.entry * has no rule here (annotate mode or a split '
                                     'child: the entry set is primary)')
                entry = rule
            else:
                entry = decode_setexpr(f[4], derive.entry_base(segs, seg), len(g.labels))
            seg.entry = rd._ids(entry, len(g.labels))
            seg.entry_total = len(seg.entry) if f[5] == '*' else parse_int(f[5])
            if f[10] == '*':
                seg.walk = None
            else:
                seg.walk = '' if f[10] == '.' else _fixed(f[10], 'bases')
                if len(seg.walk) != length:
                    raise CodecError('%d bases for length_bp %d' % (len(seg.walk), length))
            has = seg.walk is not None
            if sequences is None:
                sequences = has
            elif sequences != has:
                raise CodecError('bases on some segments only (output.sequences is one '
                                 'value per retrieval)')
            if f[9] != '*':
                if has or len(f[9]) != 1 or length == 0:
                    raise CodecError('first_base %s is written only when the bases are '
                                     'absent' % _tok(f[9]))
                seg.first_base = f[9]
            elif not has and len(parents) == 1:
                # a split child: its base is splits[].branches[].char
                raise CodecError('a split child without bases carries its first base')
        except CodecError as e:
            raise GraphletFormatError(line, str(e)) from None
        # P: base = the previous P, the first P against the entry set
        base = seg.entry
        last_to = seg.from_bp
        for pline, p in sr['P']:
            try:
                pf, pt = parse_int(p[1]), parse_int(p[2])
                if pf != last_to or pt <= pf:
                    raise CodecError('P [%d, %d) does not continue the segment at %d'
                                     % (pf, pt, last_to))
                labels = rd._ids(decode_setexpr(p[4], list(base), len(g.labels)),
                                 len(g.labels))
                run = PresenceRun(pf, pt, parse_int(p[3]), labels)
            except CodecError as e:
                raise GraphletFormatError(pline, str(e)) from None
            seg.presence.append(run)
            base = labels
            last_to = pt
        if seg.presence and last_to != seg.end_bp:
            raise GraphletFormatError(line, 'P runs end at %d, the segment at %d'
                                      % (last_to, seg.end_bp))
        seg.events = sr['E']
        for ev in seg.events:
            if ev.type == 'revisit' and ev.segment >= n:
                raise GraphletFormatError(line, 'revisit of segment %d' % ev.segment)
        try:
            if f[6] == '*':
                end = derive.rule_end(g.mode, seg, seg.entry, seg.presence, anchored[sid],
                                      arm.runs)
            else:
                end = decode_setexpr(f[6], list(seg.entry), len(g.labels))
            seg.end = rd._ids(end, len(g.labels))
        except CodecError as e:
            raise GraphletFormatError(line, str(e)) from None
        if sr['T'] is not None:
            tline, t = sr['T']
            try:
                seg.leaf = Leaf(_opt(t[1], _reason_code), _parse_extras(t[2], len(g.labels)))
            except CodecError as e:
                raise GraphletFormatError(tline, str(e)) from None
        if sr['C'] is not None:
            cline, c = sr['C']
            try:
                if seg.leaf.path_reason is None or not (
                        seg.leaf.path_reason == 'X' or seg.leaf.path_reason in RESOURCE_CODES):
                    raise CodecError('a continuation belongs to a leaf ended by the radius or '
                                     'a resource reason')
                n_bp = parse_int(c[1])
                seq = None
                if len(c) == 6:
                    # explicit when G carries no bases, and wherever the derived
                    # spelling would not reproduce the walker's (the writer's fallback)
                    seq = _fixed(c[5], 'continuation')
                    if len(seq) != n_bp or n_bp == 0:
                        raise CodecError('C sequence of %d bases for n %d' % (len(seq), n_bp))
                elif seg.walk is None and n_bp > 0:
                    raise CodecError('C carries its sequence when G carries no bases')
                seg.leaf.continuation = LeafContinuation(
                    n_bp, parse_num(c[2]), parse_int(c[3]),
                    rd._ids(decode_setexpr(c[4], list(seg.end), len(g.labels)), len(g.labels)), seq)
            except CodecError as e:
                raise GraphletFormatError(cline, str(e)) from None
        segs.append(seg)
    # derived structure
    for s in segs:
        for p in s.parents:
            segs[p].children.append(s.id)
        s.depth = segs[s.parents[0]].depth + 1 if s.parents else 0
    for s in segs:
        singles = [c for c in s.children if len(segs[c].parents) == 1]
        if singles and len(singles) < 2:
            raise GraphletFormatError(raw.seg_raw[s.id]['line'],
                                      'a split has at least two children')
        if bool(singles) != (s.split is not None):
            raise GraphletFormatError(raw.seg_raw[s.id]['line'],
                                      'G.split is 0|1 exactly on a segment that ends in a split')
        # a split head's item is gone; only its children can meet another route
        if singles and any(len(segs[x].parents) > 1 for x in s.children):
            raise GraphletFormatError(raw.seg_raw[s.id]['line'],
                                      'a segment ends in a split or a merge, not both')
    for r in arm.runs:
        if r.to_bp > segs[r.segment].end_bp or r.to_bp < segs[r.segment].from_bp:
            raise GraphletFormatError(raw.run_lines[r.id], '%s arm: run %d ends at %d '
                                      'outside its anchor %d [%d, %d]'
                                      % (arm.side, r.id, r.to_bp, r.segment,
                                         segs[r.segment].from_bp, segs[r.segment].end_bp))
    if g.mode == 'constrain':
        for leaf in derive.leaves(arm):
            seen = set()
            for e in derive.end_labels(arm, leaf):
                if e.label in seen:
                    raise GraphletFormatError(raw.seg_raw[leaf]['line'],
                                              'label %d ends twice at the leaf' % e.label)
                seen.add(e.label)
            if not seen.issuperset(segs[leaf].leaf.extras):
                raise GraphletFormatError(raw.seg_raw[leaf]['line'],
                                          'T extras for a label not alive at the leaf')
    elif arm.runs:
        raise GraphletFormatError(raw.run_lines[0], '%s arm: R records in annotate mode'
                                  % arm.side)


def _check_arm_invariants(g, arm, line):
    """What the records state about each other beyond the A counts, so that a lost
    optional record (a bin, a branch event, a continuation) is noticed too: the stored
    branch events are all of them unless the evidence boundary says otherwise (§7.2),
    the bins sum to the step counter and to the run ends (spec T25), and every leaf
    ended by the radius or a resource reason carries its continuation (finish_path)."""
    n_v = len(arm.branch_events)
    if arm.evidence_complete_to_bp is None and n_v != arm.branch_events_total:
        raise GraphletFormatError(line, '%d V records of %d branch events with no evidence '
                                  'boundary' % (n_v, arm.branch_events_total))
    if arm.evidence_complete_to_bp is not None and n_v >= arm.branch_events_total:
        raise GraphletFormatError(line, 'an evidence boundary but all %d branch events stored'
                                  % n_v)
    # no "at least one B" rule: the grammar is A B* V*, and an arm can stop before its
    # first level (MAX_STEPS tripped on the left arm censors the right root before it
    # ran a level, and in annotate mode that ends no run): its growth is [] and every sum
    # below is zero
    steps = arm.counters.get('steps')
    if steps is not None and sum(b.steps for b in arm.growth) != steps:
        raise GraphletFormatError(line, 'the B bins sum to %d steps, the counter says %d'
                                  % (sum(b.steps for b in arm.growth), steps))
    # spec T25: the bins count every split (ambiguous or a divergence) and merge
    total = {'splits': 0, 'ambiguous': 0, 'merges': 0}
    for b in arm.growth:
        total['splits'] += b.splits
        total['ambiguous'] += b.ambiguous_branches
        total['merges'] += b.reconvergences
    got = {'splits': arm.counts.splits, 'merges': arm.counts.merges,
           'ambiguous': sum(1 for s in arm.segments if s.split)}
    if total != got:
        raise GraphletFormatError(line, 'the B bins count %s, the body has %s' % (total, got))
    if g.mode == 'constrain' and arm.growth:
        ends = sum(n for b in arm.growth for n in b.label_ends.values())
        runs = sum(1 for r in arm.runs if not r.merged)
        if ends != runs:
            raise GraphletFormatError(line, 'the B bins count %d label ends, the R records %d'
                                      % (ends, runs))
    blocked = sum(1 for s in arm.segments for e in s.events if e.type == 'blocked')
    if blocked != sum(b.blocked_repeat for b in arm.growth):
        raise GraphletFormatError(line, '%d blocked events, the B bins count %d'
                                  % (blocked, sum(b.blocked_repeat for b in arm.growth)))
    for s in arm.segments:
        if g.mode == 'annotate' and s.length_bp and not s.presence:
            # every node entered is recorded (record_present): runs cover the segment
            raise GraphletFormatError(line, 'annotate segment %d has no P runs' % s.id)
        if s.leaf is not None:
            code = s.leaf.path_reason
            wants = code is not None and (code == 'X' or code in RESOURCE_CODES)
            if wants != (s.leaf.continuation is not None):
                raise GraphletFormatError(line, 'leaf %d %s a continuation' % (
                    s.id, 'lacks' if wants else 'has'))


def _parse_extras(tok, n_labels):
    out = {}
    if tok == '.':
        return out
    for item in tok.split(','):
        f = item.split(':')
        if len(f) != 4:
            raise CodecError('T extra %s' % _tok(item))
        label = parse_int(f[0])
        if label >= n_labels or label in out:
            raise CodecError('T extra label %d' % label)
        out[label] = (parse_num(f[1]), parse_int(f[2]), parse_int(f[3]))
    return out


# Top-level names of a J object that are the library's own, never a response's: a view's
# selector and the retrieval a view or a derived graphlet came from (_attach_j reads them
# back). A response carrying one would be misread on load, so it is refused
# (from_response, standalone_text).
RESERVED_ENVELOPE_NAMES = ('view', 'derived_from')


def seed_envelope(envelope, seed_index):
    """The envelope a single seed's graphlet keeps (a saved file's J line, a store entry):
    the response envelope with its `usage` reduced to the totals plus THIS seed's per_seed
    entry (the one whose index is |seed_index|, the H record's seed_index) — the other seeds'
    entries describe graphlets this one is not. Idempotent; an envelope without a per_seed
    list is returned as it is (a copy)."""
    env = dict(envelope or {})
    usage = env.get('usage')
    if isinstance(usage, dict) and isinstance(usage.get('per_seed'), list):
        usage = dict(usage)
        usage['per_seed'] = [e for e in usage['per_seed']
                             if isinstance(e, dict) and e.get('index') == seed_index]
        env['usage'] = usage
    return env


def _refuse_reserved(response):
    for name in RESERVED_ENVELOPE_NAMES:
        if name in response:
            raise ValueError('the response has a top-level %r, a name the library reserves '
                             'for its own envelope (a view or a derived graphlet): refused '
                             'rather than misread' % name)


def _attach_j(g, j):
    # 'view' and 'derived_from' are the library's own names (RESERVED_ENVELOPE_NAMES): a
    # response that carried one was refused before it could be saved
    results = j.get('results')
    if not isinstance(results, list) or len(results) != 1 or not isinstance(results[0], dict):
        raise GraphletFormatError(2, 'J: results[] reduced to this seed\'s summary')
    env = {k: v for k, v in j.items() if k not in ('results', 'view', 'derived_from')}
    # a body-only view or derived graphlet is saved with an empty envelope and summary:
    # read back as none. Any other J line keeps what it holds, an empty side included --
    # read as none, it would lose the other side and dump() would drop the J line: the
    # round trip would not be stable
    body_only = not results[0] and not env and (j.get('view') is not None
                                                or j.get('derived_from') is not None)
    g.seed_summary = None if body_only else results[0]
    g.envelope = None if body_only else env
    g.view = j.get('view')
    g.derived_from = j.get('derived_from')


def utf8_bytes(text):
    """|text| (a str body) as UTF-8 bytes. A str can hold what no UTF-8 document can: a
    lone surrogate (JSON's '\\udc80' decodes to one) -> GraphletFormatError naming its
    line, never a bare UnicodeEncodeError (a ValueError the tool layer would report as
    the agent's bad argument). The one byte count of a body: parse, from_response, the
    store and traverse_fetch use it."""
    if text.isascii():
        return text.encode('ascii')
    try:
        return text.encode('utf-8')
    except UnicodeEncodeError as e:
        raise GraphletFormatError(text.count('\n', 0, e.start) + 1,
                                  'not UTF-8: a lone surrogate (%s)'
                                  % ascii(text[e.start:e.end])) from None


def parse(text, *, budget=None):
    """-> Graphlet (BODY ONLY unless the text carries a J line: no seed_id, annotation
    counters or replay strategy; to_json(), next_request() and summary() raise
    MissingEnvelope then). Raises GraphletFormatError(line_no, msg).

    budget= (local limits): a body whose line count alone exceeds the work budget, or whose
    line list alone exceeds the memory budget, is refused before any record is read; a
    parse that stops raises LocalBudgetExceeded (not GraphletFormatError): no model, the
    text untouched -- a cut-off parse never makes a shallower graphlet (§14)."""
    b = _B.resolve(budget)
    if b is None:
        return _parse(text, None)
    with b.scope('parse', ('export_mgt',)):
        try:
            return _parse(text, b)
        except LocalBudgetExceeded as e:
            raise e.restate(of=b.done.get('of'))


def _parse(text, b):
    if b is not None:
        # admission, before anything is read: every line is parsed and the line list is
        # held, so a body whose line count alone exceeds the work budget, or whose line
        # list alone exceeds the account, is refused at once (certain lower bounds)
        n = text.count(b'\n' if isinstance(text, (bytes, bytearray)) else '\n')
        b.done = {'lines': 0, 'of': n}
        b.phase = 'admission'
        b.require(W_LINE * n, len(text) + (STR + LIST_ITEM) * n + LIST)
        b.phase = 'setup'
    if isinstance(text, (bytes, bytearray)):
        if b is not None:
            b.charge(len(text) >> 10, STR + len(text))
        try:
            text = bytes(text).decode('utf-8')
        except UnicodeDecodeError as e:
            raise GraphletFormatError(0, 'not UTF-8: %s' % e) from None
    else:
        if b is not None:
            b.charge(len(text) >> 10, len(text) + 33)
            b.release(len(text) + 33)
        utf8_bytes(text)
    if not text.endswith('\n'):
        raise GraphletFormatError(0, 'truncated document: it does not end with a line feed')
    if b is not None:
        # the digest below: C-speed hashing and one chunk encoded at a time
        chunk = 2 * (STR + min(len(text), _DIGEST_CHUNK))
        b.charge(len(text) >> 7, chunk)
        b.release(chunk)
    digest = body_digest(text)
    if b is not None:
        b.phase = 'lines'
        # the modelled account of the model as it is built (the line list included),
        # spread over the line blocks; the label sets interned are charged as they are
        model = _MODEL_PER_BYTE * len(text) + _MODEL_PER_LINE * n
        for tag, size in _MODEL_PER_RECORD:
            model += size * text.count(tag)
        b.charge(len(text) >> 10)
    # LF is the only record separator: names may hold VT, FF, NEL, LS, PS raw, so
    # str.splitlines() would cut records
    lines = text.split('\n')
    lines.pop()
    rd = _Reader(lines, b)
    if b is not None:
        rd.per_line_work = (len(text) >> 3) // max(1, n)
        rd.per_line_bytes = model / max(1, n)
    g = rd.run()
    g.body_digest = digest
    if g.has_j:
        # a saved file's coordinates (in its J line's summary) are checked against the
        # body now that its R records are read (eagerly, never at first use)
        coords.attach(g, b)
    return g


_DIGEST_CHUNK = 1 << 20


def body_digest(text):
    """The sha-256 (hex) of an MGT text's records without its J line and its Z record --
    Graphlet.body_digest, what binds a resume token (local limits): the J line holds the
    envelope and a saved view, which no list reads, and the Z record counts the J line,
    so the server's body and the standalone file saved from it have one digest. Hashed a
    chunk at a time, so a large body is never copied whole."""
    first = text.find('\n') + 1
    start = first
    if text.startswith('J ', first):
        start = text.find('\n', first) + 1
    end = text.rfind('\n', 0, len(text) - 1) + 1
    h = hashlib.sha256(text[:first].encode('utf-8'))
    for i in range(start, end, _DIGEST_CHUNK):
        h.update(text[i:min(end, i + _DIGEST_CHUNK)].encode('utf-8'))
    return h.hexdigest()


# ======================================================================= writing

def _k_or_star(v):
    return '*' if v is None else encode_kvalue(v)


def _set(ids):
    return list(ids)


def dump(g, envelope=None, *, budget=None):
    """Canonical text (§2.3). |envelope|: None -> a J line iff the graphlet was read
    with one; True -> whenever an envelope is attached; False -> body only. budget=
    (local limits): every record is charged; a stop raises LocalBudgetExceeded and no text
    is returned (done: the records written)."""
    b = _B.resolve(budget)
    if b is None:
        return _dump(g, envelope, None)
    with b.scope('dump', ('export_mgt',)):
        done = [0, 0]
        try:
            return _dump(g, envelope, b, done)
        except LocalBudgetExceeded as e:
            raise e.restate(records=done[0])


def _dump(g, envelope, b, done=None):
    out = []
    if b is None:
        w = out.append
    else:
        b.phase = 'records'

        pending = [0]

        def w(line):
            # records are charged as soon as they are built (a text's size is known only
            # then), _LINE_BLOCK at a time: the account may run up to that many records
            # ahead of the check
            n = len(line)
            done[1] += n + 1
            pending[0] += TEXT_LINE + n
            out.append(line)
            if len(out) % _LINE_BLOCK == 0:
                done[0] = len(out)
                b.charge(0, pending[0])
                pending[0] = 0
        b.charge(W_DUMP_LINE * (4 + len(g.labels) + len(g.dropped) + len(g.limitations)))
    w(' '.join([
        'H', 'mgt', str(FORMAT_VERSION), fmt_int(g.k), g.regime, g.alphabet,
        _inv(_MODE)[g.mode], _inv(_SUPPORT)[g.support], _inv(_RECONVERGE)[g.reconverge],
        fmt_int(g.cap), fmt_int(g.continuation_bp), fmt_int(g.seed_index), 'walk',
        g.index_ns or '*', g.index_fp or '*', g.index_meta_fp or '*']))
    with_j = g.has_j if envelope is None else bool(envelope)
    if with_j and (g.has_envelope or g.view is not None or g.derived_from is not None):
        w('J ' + j_text(g))
    s = g.seed
    w('S %s %d %d %d %s' % (s.validated_seed_id or '*', s.length_bp, s.num_kmers,
                            s.num_seed_labels, s.sequence))
    for d in g.dropped:
        w('X %s %s %s' % (d.reason, encode_pairs(d.runs), pct_escape(d.name)))
    prev = ''
    for lab in g.labels:
        w('L %s %d %s %s' % (_inv(_KIND)[lab.kind], lab.column,
                             '*' if lab.seq_id is None else fmt_int(lab.seq_id),
                             front_encode(prev, lab.name)))
        prev = lab.name
    o = g.outcome
    w('O %s %s %s %s' % (_inv(_O_WALKS)[o.walks], _inv(_O_BRANCH)[o.branch_diagnostics],
                         _inv(_O_LABEL)[o.label_evidence], _inv(_O_DELIVERY)[o.delivery]))
    q = g.resource_stop
    if q is not None:
        w(' '.join(['Q', q.scope, q.resource, q.phase] +
                   [_k_or_star(v) for v in (q.requested, q.effective, q.used, q.remaining)] +
                   [','.join(q.actions) if q.actions else '.', pct_escape(q.message)]))
    for lim in g.limitations:
        extra = ','.join('%s=%s' % (n, encode_kvalue(v)) for n, v in lim.extra) or '.'
        w(' '.join(['K', {'left': 'l', 'right': 'r', None: '*'}[lim.arm], lim.kind, lim.knob,
                    encode_kvalue(lim.limit), encode_kvalue(lim.observed),
                    _star(lim.complete_to_bp, fmt_int), extra, pct_escape(lim.effect)]))
    for side in ('left', 'right'):
        arm = g.arms.get(side)
        if arm is not None:
            _dump_arm(g, arm, w, b)
    n = len(out) + 1
    if b is not None:
        done[0] = len(out)
        b.charge(0, pending[0])
        b.phase = 'join'
        b.charge(done[1] >> 8, STR + done[1] + 16)
    out.append('Z %d' % n)
    # the final line feed joined with the rest: adding it to the joined text copied the
    # whole text once more at the peak
    out.append('')
    return '\n'.join(out)


def _dump_arm(g, arm, w, bud=None):
    counters = ','.join('%s=%s' % (k, fmt_int(v)) for k, v in arm.counters.items()) or '.'
    c = arm.cap_trigger
    cap = '*' if c is None else ','.join([c.reason, fmt_int(c.at_bp), fmt_int(c.segment),
                                          fmt_int(c.live_paths), fmt_int(c.live_labels),
                                          fmt_bool(c.exact), fmt_num(c.demand)])
    k = arm.counts
    w(' '.join([
        'A', arm.code, _inv(_STATUS)[arm.status], fmt_int(arm.complete_to_bp),
        _inv(_SCOPE)[arm.scope], fmt_int(arm.frontier.live_paths),
        fmt_int(arm.frontier.live_labels), fmt_bool(arm.frontier.exact),
        fmt_int(arm.labels_per_node.max_seen), fmt_int(arm.labels_per_node.nodes_truncated),
        counters, cap, fmt_int(arm.branch_events_total),
        _star(arm.evidence_complete_to_bp, fmt_int), fmt_int(k.segments), fmt_int(k.runs),
        fmt_int(k.leaves), fmt_int(k.splits), fmt_int(k.merges), fmt_int(k.bases)]))
    for b in arm.growth:
        ends = ','.join('%s:%d' % (code, b.label_ends[code]) for code in REASON_ORDER
                        if b.label_ends.get(code)) or '.'
        w(' '.join(['B'] + [fmt_int(x) for x in (b.from_bp, b.max_live_paths,
                                                 b.distinct_live_labels, b.live_pairs)] +
                   [fmt_bool(b.exact)] +
                   [fmt_int(x) for x in (b.steps, b.divergences, b.ambiguous_branches,
                                         b.splits, b.reconvergences, b.bubbles, b.tips,
                                         b.blocked_repeat)] + [ends]))
    for be in arm.branch_events:
        refused = ';'.join('%s:%s:%s' % (r.char, r.cause, encode_ranges(r.labels))
                           for r in be.refused) or '.'
        w(' '.join(['V', fmt_int(be.at_bp), fmt_int(be.segment), be.chars or '.',
                    ','.join(fmt_int(x) for x in be.labels_per_successor) or '.',
                    encode_ranges(be.ambiguous), encode_ranges(be.dropped), refused]))
    segs = arm.segments
    rbs = derive.runs_by_segment(arm)
    nsl = g.seed.num_seed_labels
    if bud is not None:
        bud.uses((id(arm), 'dump_prices'), 3 * len(arm.segments) + len(arm.runs),
                 16 * len(arm.segments))
        pw, pt = _dump_prices(arm, rbs)
        bud.charge(W_DUMP_LINE * (len(arm.growth) + len(arm.branch_events) + len(arm.runs))
                   + arm.cache['dump_branch_ids'])
    if bud is not None:
        tmp = 0
    for s in segs:
        if bud is not None and s.id % _LINE_BLOCK == 0:
            # the next block of segments' records: their sets encoded (each against its
            # base, the '*' rules recomputed), P runs, events, leaves and continuations;
            # and the temporaries of the widest set among them, dropped after the block
            bud.release(tmp)
            j = s.id + _LINE_BLOCK
            tmp = max(pt[s.id:j])
            bud.charge(sum(pw[s.id:j]), tmp)
        entry = _set(s.entry)
        rule = derive.rule_entry(g.mode, nsl, s)
        if rule is not None and rule == entry:
            entry_tok = '*'
        else:
            entry_tok = encode_setexpr(entry, derive.entry_base(segs, s))
        end = _set(s.end)
        if derive.rule_end(g.mode, s, s.entry, s.presence, rbs[s.id], arm.runs) == end:
            end_tok = '*'
        else:
            end_tok = encode_setexpr(end, entry)
        if all(len(p) == 0 for p in s.partition) and \
                len(s.partition) == (len(s.parents) if len(s.parents) > 1 else 0):
            part_tok = '*'
        else:
            part_tok = '|'.join(encode_ranges(p) for p in s.partition)
        w(' '.join([
            'G', ','.join(fmt_int(p) for p in s.parents) or '*', fmt_int(s.from_bp),
            fmt_int(s.length_bp), entry_tok,
            '*' if s.entry_total == len(s.entry) else fmt_int(s.entry_total),
            end_tok, part_tok, _star(s.split, fmt_bool),
            # only a split child needs it (its branch's char), and only without bases
            s.first_base if s.walk is None and len(s.parents) == 1 and s.first_base else '*',
            '*' if s.walk is None else (s.walk or '.')]))
        base = entry
        for p in s.presence:
            labels = _set(p.labels)
            w('P %d %d %d %s' % (p.from_bp, p.to_bp, p.total, encode_setexpr(labels, base)))
            base = labels
        for ev in s.events:
            w(_event_line(ev))
        if s.leaf is not None:
            extras = ','.join('%d:%s:%d:%d' % (l, fmt_num(x[0]), x[1], x[2])
                              for l, x in sorted(s.leaf.extras.items())
                              if x[0] or x[1] or x[2]) or '.'
            w('T %s %s' % (s.leaf.path_reason or '*', extras))
            c = s.leaf.continuation
            if c is not None:
                line = 'C %d %s %d %s' % (c.n, fmt_num(c.loss_used), c.branches_used,
                                          encode_setexpr(_set(c.labels), end))
                if c.sequence is not None and c.n > 0 and (
                        s.walk is None or c.sequence != _derived_continuation(g, arm, s)):
                    line += ' ' + c.sequence
                w(line)
    if bud is not None:
        bud.release(tmp)
    for r in arm.runs:
        fields = ['R', fmt_int(r.segment), fmt_int(r.label), fmt_int(r.from_bp),
                  fmt_int(r.to_bp), r.end, fmt_int(r.route_bp),
                  '*' if r.from_label is None else '%d:%s' % (r.from_label, fmt_num(r.cost)),
                  _star(r.prev_run, fmt_int), _star(r.structural_successors, fmt_int),
                  fmt_int(r.branches), fmt_num(r.loss)]
        if r.code == 'B':
            fields.append(fmt_num(r.needed_budget))
        w(' '.join(fields))


def _dump_prices(arm, rbs):
    """Per segment the work of writing its records and the temporaries of its widest set
    (measured: about four sets and five int lists of it alive at once while the entry and
    end sets are encoded and their '*' rules recomputed): facts of the model, cached."""
    got = arm.cache.get('dump_prices')
    if got is None:
        pw, pt = [], []
        for s in arm.segments:
            ids = 2 * (len(s.entry) + len(s.end)) + len(rbs[s.id])
            for p in s.presence:
                ids += len(p.labels)
            for p in s.partition:
                ids += len(p)
            lines = 1 + len(s.presence) + len(s.events) + (2 if s.leaf is not None else 0)
            pw.append(W_DUMP_LINE * lines + ids)
            wide = max([len(s.entry), len(s.end)] + [len(p.labels) for p in s.presence])
            pt.append(4 * (216 + 54 * wide) + 5 * 36 * wide)
        got = arm.cache['dump_prices'] = (pw, pt)
        arm.cache['dump_branch_ids'] = sum(
            len(be.ambiguous) + len(be.dropped) + sum(len(r.labels) for r in be.refused)
            for be in arm.branch_events)
    return got


def _derived_continuation(g, arm, seg):
    c = seg.leaf.continuation
    saved, c.sequence = c.sequence, None
    try:
        return derive.continuation_sequence(g, arm, seg.id)
    finally:
        c.sequence = saved


def _event_line(ev):
    t = ev.type
    if t == 'switch':
        return 'E %d s %d %d %s' % (ev.at_bp, ev.from_label, ev.to_label, fmt_num(ev.cost))
    if t == 'blocked':
        return 'E %d b %s %s %d %s' % (ev.at_bp, ev.char, ev.reason, ev.total,
                                       encode_ranges(ev.labels))
    if t == 'hairpin':
        return 'E %d h %s %d %s %s' % (ev.at_bp, ev.char, ev.total, encode_ranges(ev.labels),
                                       'f' if ev.followed else 's')
    if t == 'revisit':
        return 'E %d v %d %s' % (ev.at_bp, ev.segment,
                                 '=' if ev.delta is None else fmt_int(ev.delta))
    if t == 'tip':
        return 'E %d t %s %d' % (ev.at_bp, ev.char, ev.length_bp)
    if t == 'bubble':
        return 'E %d u %d %s' % (ev.at_bp, ev.length_bp, ev.alleles)
    raise ValueError('event type %r' % t)


def j_object(g):
    """The J line's object: the response envelope with results[] reduced to this seed's
    summary (the graphlet string removed) and its usage to the totals plus this seed's
    per_seed entry (seed_envelope), plus a view's selector when it is one."""
    j = seed_envelope(g.envelope, g.seed_index)
    summary = {k: v for k, v in (g.seed_summary or {}).items()
               if k not in ('graphlet', 'graphlet_bytes', 'graphlet_lines')}
    j['results'] = [summary]
    if g.view is not None:
        j['view'] = g.view
    if g.derived_from is not None:
        j['derived_from'] = g.derived_from
    return j


def j_text(g):
    # one line, deterministic: sorted keys, no whitespace, raw UTF-8 (json escapes every
    # control character, so the line cannot break)
    return json.dumps(j_object(g), sort_keys=True, separators=(',', ':'),
                      ensure_ascii=False)


def is_canonical(text):
    """True iff |text| is a valid document in canonical form (dump(parse(text)) == text);
    every writer output must be (T37)."""
    try:
        return dump(parse(text)) == text
    except GraphletFormatError:
        return False


# ======================================================================= envelope

def _transport_count(result, key):
    """result[key] (graphlet_bytes, graphlet_lines) as an int, None when absent. Anything
    else is a GraphletFormatError: a str, None or list would compare unequal and then fail
    the '%d' of the message with a TypeError, and True would pass as 1."""
    if key not in result:
        return None
    v = result[key]
    if isinstance(v, bool) or not isinstance(v, int):
        raise GraphletFormatError(0, '%s is not an integer: %r' % (key, v))
    return v


def _check_transport(body, result):
    """The body against its summary's byte and line counts (a body cut in transport),
    the cheap checks before any parse: the one copy of them (from_response(),
    standalone_text() and the store's unparsed put call it)."""
    want = _transport_count(result, 'graphlet_bytes')
    if want is not None:
        # an ASCII body's bytes are its characters: no encoded copy to count them
        n = len(body) if body.isascii() else len(utf8_bytes(body))
        if n != want:
            raise GraphletFormatError(0, 'the body has %d bytes, graphlet_bytes says %d '
                                      '(truncated in transport)' % (n, want))
    want = _transport_count(result, 'graphlet_lines')
    if want is not None and body.count('\n') != want:
        raise GraphletFormatError(0, 'graphlet_lines says %d, the body has %d'
                                  % (want, body.count('\n')))


def from_response(result, response, *, budget=None):
    """results[i] of a detail: graphlet response + the response (the envelope) ->
    Graphlet with the summary attached. A derivation failure ({seed, error}) has no
    graphlet and raises ValueError, as does a response with a reserved top-level name
    (RESERVED_ENVELOPE_NAMES)."""
    body = result.get('graphlet')
    if body is None:
        raise ValueError('this result carries no graphlet%s' % (
            ': ' + result['error'] if 'error' in result else ''))
    _refuse_reserved(response)
    b = _B.resolve(budget)
    if b is not None:
        with b.scope('parse', ('export_mgt',)):
            if 'graphlet_bytes' in result:
                b.charge(len(body) >> 10, 4 * len(body) + 33)
                b.release(4 * len(body) + 33)
            return _from_response_checked(result, response, body, b)
    return _from_response_checked(result, response, body, None)


def _from_response_checked(result, response, body, b):
    _check_transport(body, result)
    g = parse(body, budget=b)
    g.seed_summary = {k: v for k, v in result.items() if k != 'graphlet'}
    g.envelope = {k: v for k, v in response.items() if k != 'results'}
    # the record coordinates of the summary against the body, eagerly: an inconsistent
    # block is a GraphletFormatError here, like a cut body
    coords.attach(g, b)
    return g


def standalone_text(body, result, response):
    """The standalone .mgt text of results[i] of a detail: graphlet response — H, the J line
    (the envelope with this seed only, seed_envelope) and the body — spliced WITHOUT parsing
    the body: for a server body (canonical, T37) it is byte for byte
    dump(from_response(result, response), envelope=True), at the cost of a copy, so that a
    large retrieval can be saved or spooled without building its model. The body is checked
    against graphlet_bytes and graphlet_lines, and for its H and Z records, but NOT
    validated: load() or parse() of the text does that (a body that is not canonical gives a
    text that is not what dump() would write). Raises ValueError for a body that already
    carries a J line (not a server body) and for a response with a reserved top-level name,
    GraphletFormatError for a body cut in transport or without its H / Z record."""
    _refuse_reserved(response)
    if isinstance(body, (bytes, bytearray)):
        try:
            body = bytes(body).decode('utf-8')
        except UnicodeDecodeError as e:
            raise GraphletFormatError(0, 'not UTF-8: %s' % e) from None
    _check_transport(body, result)
    if not body.endswith('\n'):
        raise GraphletFormatError(0, 'truncated document: it does not end with a line feed')
    first = body.find('\n')
    head = body[:first].split(' ')
    if len(head) != 16 or head[0] != 'H' or head[1] != 'mgt':
        raise GraphletFormatError(1, 'expected the H record')
    try:
        seed_index = parse_int(head[11])
    except CodecError as e:
        raise GraphletFormatError(1, 'H seed_index: %s' % e) from None
    # tested in place: a slice of the rest of the body only to read its first two
    # characters would copy the whole body once more, held until the splice
    if body.startswith('J ', first + 1):
        raise ValueError('the body already carries a J line: a saved file, not a server body')
    last = body.rfind('\n', 0, len(body) - 1) + 1
    z = body[last:-1].split(' ')
    if len(z) != 2 or z[0] != 'Z':
        raise GraphletFormatError(body.count('\n'), 'expected the Z record')
    try:
        lines = parse_int(z[1])
    except CodecError as e:
        raise GraphletFormatError(body.count('\n'), 'Z: %s' % e) from None
    envelope = seed_envelope({k: v for k, v in response.items() if k != 'results'},
                             seed_index)
    envelope['results'] = [{k: v for k, v in result.items()
                            if k not in ('graphlet', 'graphlet_bytes', 'graphlet_lines')}]
    j = json.dumps(envelope, sort_keys=True, separators=(',', ':'), ensure_ascii=False)
    # the Z record counts the lines, the J line now among them
    return ''.join((body[:first + 1], 'J ', j, '\n', body[first + 1:last],
                    'Z %d\n' % (lines + 1)))


def save(g, path, *, budget=None):
    """Write the standalone .mgt file: H, J (envelope with this seed only), body.
    Atomic: a reader never sees a partial file. budget= (local limits): the text is made
    in full before the file is opened, so a stop writes no file."""
    b = _B.resolve(budget)
    if b is not None:
        with b.scope('save', ('export_mgt',)):
            text = dump(g, envelope=True, budget=b)
            b.phase = 'write'
            # the text's UTF-8 bytes, written without a buffer (_write_text)
            b.charge(len(text) >> 8, 2 * (len(text) + 33))
            return _write_text(text, path)
    text = dump(g, envelope=True)
    return _write_text(text, path)


def _write_text(text, path):
    data = text.encode('utf-8')
    d = os.path.dirname(os.path.abspath(path))
    fd, tmp = tempfile.mkstemp(prefix='.mgt-', dir=d)
    try:
        _write_all(fd, data)
        os.replace(tmp, path)
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise
    return len(data)


def _write_all(fd, data, sync=False):
    """Write the bytes |data| in full to the file descriptor |fd| and close it -- without a
    buffer: they go out in one write, and a buffered file holds io.DEFAULT_BUFFER_SIZE (128
    KiB from Python 3.14, more on a file system with larger blocks) that the local-limits
    account of a save does not charge, which would put it below the traced peak on bodies of
    10-40 KB. |sync|: the data is flushed to the device (fsync) before the file is closed --
    what a write renamed into place needs so that a crash does not leave the new name on an
    empty or partly written file (the store's spool)."""
    with os.fdopen(fd, 'wb', buffering=0) as f:
        view = memoryview(data)
        while view:
            # a raw write may write less than it was given
            n = f.write(view)
            if not n:
                raise OSError('%s accepted none of %d bytes' % (f.name, len(view)))
            view = view[n:]
        if sync:
            os.fsync(f.fileno())


def load(path, *, budget=None):
    """A saved .mgt file -> its Graphlet (or the GraphletView it saved). budget=
    (local limits): the file's text, the parse and the view are charged; a stop raises
    LocalBudgetExceeded and the file is left as it is."""
    b = _B.resolve(budget)
    if b is None:
        with open(path, 'r', encoding='utf-8', newline='') as f:
            text = f.read()
        g = parse(text)
        if g.view is not None:
            from . import ops
            return ops.view_from_saved(g)
        return g
    with b.scope('load', ('export_mgt',)):
        n = os.path.getsize(path)
        b.charge(n >> 10, STR + 4 * n)
        with open(path, 'r', encoding='utf-8', newline='') as f:
            text = f.read()
        g = parse(text, budget=b)
        del text
        if g.view is not None:
            from . import ops
            return ops.view_from_saved(g, budget=b)
        return g


def body_text(g):
    """The body as the server sends it (no J)."""
    return dump(g, envelope=False)


def labels_array(ids):
    return array('I', ids)


def code_of_qual(code, qualifier):
    return end_token(code, qualifier)


QUALIFIERS = dict(QUAL)
