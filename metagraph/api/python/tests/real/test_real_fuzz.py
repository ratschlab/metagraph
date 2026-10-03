"""Mutation fuzzing of REAL graphlet bodies (DESIGN-traverse-graphlet.md §2, §6, §14).

The bodies are the server's own MGT documents from the retrieval cache: per index and
strategy the smallest cached body (so every record type the server writes is covered:
G/P/E/T/C/R/V/B/K/X/L, constrain and annotate, merge and keep, trace, no_sequences, cut
lists, cut events, time budgets, batches). Each is mutated deterministically:

  must raise GraphletFormatError
    truncate_lines   the body cut after line j (every j on small bodies)
    truncate_bytes   cut inside a line (no final LF)
    drop_line        one record removed
    dup_line         one record duplicated
    count_flip       Z or an A count (segments, runs, leaves, splits, merges, bases) +-1
    huge_range       a label RANGES / SETEXPR set to 0-2^64-1, 0-4294967295, ...
    u64_overflow     an integer field set to 2^64 or 25 digits
    huge_digits      an integer field set to 5000 digits                 (PRODUCT BUG 1)
    bad_utf8         invalid UTF-8 bytes in a label name (bytes input)
    lone_surrogate   a lone surrogate in a label name (str input)      (PRODUCT BUG 2)
    crlf             CR LF line ends
  may parse (then dump/parse must round-trip canonically) or raise GraphletFormatError
    swap_lines       two records swapped
    field_corrupt    one field replaced by a hostile token

The parser never hangs, crashes or allocates without bound: the parse fuzz runs in a
CHILD process that the parent watches (wall clock, per-case time, RSS via ps; it kills
the child and names the case on a breach), and every case's time and tracemalloc peak
is bounded by a multiple of the unmutated body's own parse.

The tool layer (GraphletTools / GraphletStore) gets the same corrupt bodies through
traverse_fetch (a stub transport returning the real response with the mutated body),
graphlet_load (mutated saved files) and GraphletStore.put: a structured format_error,
nothing stored, no content of the file echoed.

    cd api/python/tests/real && python3 -m unittest -v test_real_fuzz
"""

import copy
import hashlib
import json
import os
import random
import re
import resource
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import tracemalloc
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import realdata as R  # noqa: E402  (puts api/python on sys.path)

from metagraph.traverse import (  # noqa: E402
    GraphletFormatError, GraphletStore, TraverseClient, dump, parse,
)
from metagraph.traverse.mcp_tools import GraphletTools  # noqa: E402

FUZZ_MAX_BODY = 64 << 10       # bodies fuzzed: per index x strategy the smallest (of at
FUZZ_MIN_BODY = 1000           # least 1 KB) and the largest (of at most 64 KB)
WALL_CAP_S = 900.0             # the whole child run
CASE_CAP_S = 20.0              # one case (a hang)
RSS_CAP_BYTES = 1536 << 20     # the child's resident set (runaway allocation)
CASE_TIME = (0.5, 25.0)        # a case may take max(0.5 s, 25 x the unmutated parse)
CASE_PEAK = (8 << 20, 4.0)     # ... and allocate max(8 MB, 4 x the unmutated parse's peak)
CASE_RSS_GROWTH = 64 << 20     # growth of the child's max RSS within one case

MUST = ('truncate_lines', 'truncate_bytes', 'drop_line', 'dup_line', 'count_flip',
        'huge_range', 'u64_overflow', 'huge_digits', 'bad_utf8', 'lone_surrogate', 'crlf')
MAY = ('swap_lines', 'field_corrupt')
FAMILIES = MUST + MAY
# the families that fail today, each with the product bug it demonstrates. Fixed
# (2026-10-02), both now must raise GraphletFormatError naming a line:
#   huge_digits     BUG-1: an integer token of more than 4300 digits raised a bare
#                   ValueError (Python's int-string conversion limit) from parse_int /
#                   decode_ranges; the integer tokens are bounded to 20 digits now, and
#                   an error message repeats at most 32 characters of a token
#   lone_surrogate  BUG-2: a lone surrogate in a str body (what JSON "\\udc80" decodes
#                   to) raised UnicodeEncodeError from front_decode_parts (parse) and from
#                   the graphlet_bytes check (from_response)
PRODUCT_BUGS = {}

# the R-record range checks of parser._resolve_arm (BUG-5, fixed: they name the R
# record's line and the arm)
RUN_RANGE_ERROR = re.compile(r'(?:line \d+: )?(?:(?:left|right) arm: )?run \d+ ends ')
# how long a format error's message may be: a token is echoed cut to 32 characters
MAX_MESSAGE = 400

HOSTILE = ['', '*', '.', '-1', 'x', '0', '1', '18446744073709551615', '18446744073709551616',
           '1.5', '01', 'nan', 'inf', '-0', '1e9', '!', '!+', '!.', '0-0', '5-3', '0,0',
           '3,1', '%ZZ', '%0', 'é', '\t', '\x00', 'A' * 300, 'i:', 's:%', 'u:', '::']
HUGE_RANGES = ['0-18446744073709551615', '0-4294967295', '0-99999999', '5-300000000',
               '!+0-18446744073709551615', '!0-99999999', '0,2-18446744073709551615']
BAD_UTF8 = [b'\xff', b'\xc3', b'\xed\xa0\x80', b'\xc0\xaf', b'\xf8\x88\x80\x80\x80',
            b'\xe2\x82']


# ------------------------------------------------------------------ the bodies

_BODIES = []


def fuzz_bodies():
    """[(index, cell, result_i, text)], deterministic: per index and strategy the
    smallest cached body of at least FUZZ_MIN_BODY bytes (else the smallest) and the
    largest of at most FUZZ_MAX_BODY; every result of a batch."""
    if _BODIES:
        return list(_BODIES)
    out = []
    size = lambda r: (r.get('graphlet') or {}).get('body_bytes') or 0
    for index in R.INDEXES:
        for st in R.STRATEGIES:
            rows = R.cells(index, strategy=st, max_body_bytes=FUZZ_MAX_BODY)
            if not rows:
                continue
            big = [r for r in rows if size(r) >= FUZZ_MIN_BODY]
            picks = [min(big or rows, key=lambda r: (size(r), r['cell'])),
                     max(rows, key=lambda r: (size(r), r['cell']))]
            for row in picks[:1] if picks[0] is picks[1] else picks:
                c = R.load_cell(index, row)
                for i in range(c.n_results()):
                    t = c.text(i)
                    if t is not None:
                        out.append((index, row['cell'], i, t))
    _BODIES.extend(out)
    return out


def _rng(*key):
    return random.Random(int(hashlib.sha256('/'.join(map(str, key)).encode())
                             .hexdigest()[:16], 16))


def _lines(text):
    lines = text.split('\n')
    assert lines[-1] == ''
    return lines[:-1]


def _join(lines):
    return '\n'.join(lines) + '\n'


def _set_field(line, idx, tok, maxsplit=None):
    f = line.split(' ') if maxsplit is None else line.split(' ', maxsplit)
    f[idx] = tok
    return ' '.join(f)


def _int_fields(lines):
    """(line, field) of integer fields in records of every kind."""
    spots = []
    for j, l in enumerate(lines):
        t = l[:2]
        if t == 'Z ':
            spots.append((j, 1))
        elif t == 'G ':
            spots += [(j, 2), (j, 3)]
        elif t == 'P ':
            spots += [(j, 1), (j, 2), (j, 3)]
        elif t == 'A ':
            spots += [(j, 3), (j, 14)]
        elif t == 'R ':
            spots += [(j, 1), (j, 3), (j, 4)]
        elif t == 'B ':
            spots += [(j, 1)]
        elif t == 'S ':
            spots += [(j, 2), (j, 3)]
    return spots


def _range_fields(lines):
    """(line, field) of label sets: G entry, P labels, E b labels, V ranges."""
    spots = []
    for j, l in enumerate(lines):
        f = l.split(' ')
        if f[0] == 'G' and len(f) == 11:
            spots.append((j, 4))
        elif f[0] == 'P' and len(f) == 5:
            spots.append((j, 4))
        elif f[0] == 'E' and len(f) == 7 and f[2] == 'b':
            spots.append((j, 6))
        elif f[0] == 'V' and len(f) == 8:
            spots += [(j, 5), (j, 6)]
    return spots


def mutations(text, family, key):
    """-> [(payload, desc)] of |family| for one body (deterministic in |key|)."""
    lines = _lines(text)
    n = len(lines)
    rng = _rng(key, family)
    out = []
    if family == 'truncate_lines':
        js = list(range(n)) if n <= 40 else sorted(set([0, 1, 2, n - 1, n - 2] + rng.sample(
            range(n), 30)))
        for j in js:
            out.append((_join(lines[:j]) if j else '', 'first %d of %d lines' % (j, n)))
    elif family == 'truncate_bytes':
        for _ in range(4):
            off = rng.randrange(1, len(text) - 1)
            while text[off - 1] == '\n':
                off -= 1
            out.append((text[:off], 'cut at char %d of %d' % (off, len(text))))
    elif family == 'drop_line':
        for j in sorted(set([0, n - 1] + rng.sample(range(n), min(n, 5)))):
            out.append((_join(lines[:j] + lines[j + 1:]), 'line %d (%s) dropped'
                        % (j + 1, lines[j][:1])))
    elif family == 'dup_line':
        for j in sorted(set([n - 1] + rng.sample(range(n), min(n, 4)))):
            out.append((_join(lines[:j + 1] + lines[j:]), 'line %d (%s) duplicated'
                        % (j + 1, lines[j][:1])))
    elif family == 'swap_lines':
        for _ in range(8):
            a, b = rng.sample(range(n), 2)
            if lines[a] == lines[b]:
                continue
            m = list(lines)
            m[a], m[b] = m[b], m[a]
            out.append((_join(m), 'lines %d (%s) and %d (%s) swapped'
                        % (a + 1, lines[a][:1], b + 1, lines[b][:1])))
    elif family == 'count_flip':
        for j, l in enumerate(lines):
            if l.startswith('Z '):
                v = int(l.split(' ')[1])
                for d in (1, -1):
                    out.append((_join(lines[:j] + ['Z %d' % (v + d)] + lines[j + 1:]),
                                'Z %+d' % d))
            elif l.startswith('A '):
                f = l.split(' ')
                for k in range(14, 20):
                    v = int(f[k])
                    d = 1 if v == 0 or rng.random() < 0.5 else -1
                    m = list(f)
                    m[k] = str(v + d)
                    out.append((_join(lines[:j] + [' '.join(m)] + lines[j + 1:]),
                                'A line %d count field %d %+d' % (j + 1, k, d)))
    elif family == 'field_corrupt':
        for _ in range(24):
            j = rng.randrange(n)
            f = lines[j].split(' ')
            if len(f) < 2:
                continue
            k = rng.randrange(1, len(f))
            tok = rng.choice(HOSTILE)
            out.append((_join(lines[:j] + [_set_field(lines[j], k, tok)] + lines[j + 1:]),
                        'line %d (%s) field %d := %r' % (j + 1, f[0], k, tok[:20])))
    elif family == 'huge_range':
        spots = _range_fields(lines)
        for j, k in rng.sample(spots, min(len(spots), 6)):
            tok = rng.choice(HUGE_RANGES)
            out.append((_join(lines[:j] + [_set_field(lines[j], k, tok)] + lines[j + 1:]),
                        'line %d (%s) field %d := %s' % (j + 1, lines[j][:1], k, tok)))
    elif family in ('u64_overflow', 'huge_digits'):
        spots = _int_fields(lines)
        toks = ['18446744073709551616', '9' * 25] if family == 'u64_overflow' else ['9' * 5000]
        for j, k in rng.sample(spots, min(len(spots), 5)):
            tok = rng.choice(toks)
            out.append((_join(lines[:j] + [_set_field(lines[j], k, tok)] + lines[j + 1:]),
                        'line %d (%s) field %d := %d digits' % (j + 1, lines[j][:1], k,
                                                                len(tok))))
        if family == 'huge_digits':
            ranges = _range_fields(lines)
            for j, k in ranges[:1]:
                out.append((_join(lines[:j] + [_set_field(lines[j], k, '0,' + '9' * 5000)]
                                  + lines[j + 1:]),
                            'line %d (%s) RANGES field %d := 0,<5000 digits>'
                            % (j + 1, lines[j][:1], k)))
    elif family == 'bad_utf8':
        data = text.encode('utf-8')
        ls = [j for j, l in enumerate(lines) if l.startswith('L ')]
        if ls:
            raw = data.split(b'\n')
            for bad in BAD_UTF8[:5]:
                j = rng.choice(ls)
                m = list(raw)
                m[j] = m[j] + bad if rng.random() < 0.5 else m[j][:-1] + bad + m[j][-1:]
                out.append((b'\n'.join(m), 'L line %d name with %r' % (j + 1, bad)))
    elif family == 'lone_surrogate':
        ls = [j for j, l in enumerate(lines) if l.startswith('L ')]
        for s in ('\udc80', '\ud800'):
            if ls:
                j = rng.choice(ls)
                out.append((_join(lines[:j] + [lines[j] + s] + lines[j + 1:]),
                            'L line %d name + %r' % (j + 1, s)))
    elif family == 'crlf':
        out.append((text.replace('\n', '\r\n'), 'CR LF line ends'))
    else:
        raise KeyError(family)
    return out


def all_cases(families=FAMILIES):
    """[(case_id, index, cell, i, family, desc, payload)] in a fixed order."""
    cid = 0
    for index, cell, i, text in fuzz_bodies():
        for family in families:
            for payload, desc in mutations(text, family, '%s/%s/%d' % (index, cell, i)):
                yield cid, index, cell, i, family, desc, payload
                cid += 1


# ------------------------------------------------------------------ the child

def _maxrss():
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r if sys.platform == 'darwin' else r * 1024


def _measure(payload):
    tracemalloc.reset_peak()
    m0 = tracemalloc.get_traced_memory()[0]
    t0 = time.perf_counter()
    try:
        g = parse(payload)
        outcome, exc, msg = 'parsed', None, None
    except GraphletFormatError as e:
        g, outcome, exc, msg = None, 'format_error', 'GraphletFormatError', str(e)[:200]
    except Exception as e:                       # noqa: BLE001 -- what the test reports
        g, outcome, exc, msg = None, 'other', type(e).__name__, str(e)[:200]
    dt = time.perf_counter() - t0
    peak = tracemalloc.get_traced_memory()[1] - m0
    return g, outcome, exc, msg, dt, peak


def child_main(out_path):
    """Run every case; one JSON line per start and per result (flushed), so that the
    parent can name a case that hangs or blows up."""
    tracemalloc.start()
    base = {}
    with open(out_path, 'w') as out:
        n = 0
        for cid, index, cell, i, family, desc, payload in all_cases():
            key = (cell, i)
            if key not in base:
                text = R.load_cell(index, cell).text(i)
                ms = [_measure(text) for _ in range(3)]
                if ms[0][1] != 'parsed':
                    raise SystemExit('the unmutated body of %s/%d does not parse' % key)
                base[key] = (min(m[4] for m in ms), max(m[5] for m in ms))
            out.write(json.dumps({'start': cid, 'cell': cell, 'i': i, 'family': family,
                                  'desc': desc}) + '\n')
            out.flush()
            rss0 = _maxrss()
            g, outcome, exc, msg, dt, peak = _measure(payload)
            roundtrip = None
            if g is not None:
                try:
                    d1 = dump(g)
                    roundtrip = dump(parse(d1)) == d1
                except Exception as e:           # noqa: BLE001
                    roundtrip = '%s: %s' % (type(e).__name__, str(e)[:120])
            del g
            out.write(json.dumps({
                'id': cid, 'index': index, 'cell': cell, 'i': i, 'family': family,
                'desc': desc, 'outcome': outcome, 'exc': exc, 'msg': msg, 'time_s': dt,
                'peak': peak, 'base_time_s': base[key][0], 'base_peak': base[key][1],
                'rss_growth': _maxrss() - rss0, 'roundtrip': roundtrip,
                'payload_bytes': len(payload)}) + '\n')
            out.flush()
            n += 1
        out.write(json.dumps({'done': n, 'max_rss': _maxrss()}) + '\n')


def _rss_of(pid):
    try:
        out = subprocess.run(['ps', '-o', 'rss=', '-p', str(pid)], capture_output=True,
                             text=True, timeout=10).stdout.strip()
        return int(out) * 1024 if out else 0
    except (OSError, ValueError, subprocess.SubprocessError):
        return 0


def run_child(workdir):
    """Run child_main in a watched subprocess -> (results, report). report: returncode,
    killed (why), the case in progress when it was killed, done line, wall time."""
    out_path = os.path.join(workdir, 'fuzz.jsonl')
    err_path = os.path.join(workdir, 'fuzz.stderr')
    env = dict(os.environ, METAGRAPH_REAL_SCRATCH=R.SCRATCH)
    with open(err_path, 'w') as err:
        proc = subprocess.Popen([sys.executable, os.path.abspath(__file__), '--fuzz-child',
                                 out_path], stdout=subprocess.DEVNULL, stderr=err, env=env,
                                cwd=HERE)
    t0 = time.time()
    killed, current, current_since, max_rss = None, None, None, 0
    pos, results, done, last_rss_check = 0, {}, None, 0.0
    partial = ''
    while True:
        rc = proc.poll()
        if os.path.exists(out_path):
            with open(out_path) as f:
                f.seek(pos)
                chunk = f.read()
                pos = f.tell()
            partial += chunk
            *complete, partial = partial.split('\n')
            for line in complete:
                if not line:
                    continue
                j = json.loads(line)
                if 'start' in j:
                    current, current_since = j, time.time()
                elif 'done' in j:
                    done = j
                else:
                    results[j['id']] = j
                    current = None
        if rc is not None:
            break
        now = time.time()
        if now - last_rss_check > 0.5:
            last_rss_check = now
            rss = _rss_of(proc.pid)
            max_rss = max(max_rss, rss)
            if rss > RSS_CAP_BYTES:
                killed = 'rss %d MB > %d MB' % (rss >> 20, RSS_CAP_BYTES >> 20)
        if now - t0 > WALL_CAP_S:
            killed = 'wall clock %.0f s > %.0f s' % (now - t0, WALL_CAP_S)
        if current is not None and current_since and now - current_since > CASE_CAP_S:
            killed = 'case ran %.0f s > %.0f s' % (now - current_since, CASE_CAP_S)
        if killed:
            proc.kill()
            proc.wait()
            break
        time.sleep(0.1)
    with open(err_path) as f:
        stderr = f.read()[-4000:]
    return results, {'returncode': proc.returncode, 'killed': killed, 'in_progress': current,
                     'done': done, 'wall_s': time.time() - t0, 'max_rss_seen': max_rss,
                     'stderr': stderr}


# ------------------------------------------------------------------ the parse fuzz

def _fmt(r):
    return '%s/%d %s [%s]: %s %s %s' % (r['cell'], r['i'], r['family'], r['desc'], r['outcome'],
                                       r.get('exc') or '', (r.get('msg') or '')[:120])


@R.skip_unless_cache()
class TestParserFuzz(unittest.TestCase):
    """The parse fuzz, run once in a watched child process; the test methods assert
    over its results."""

    results = None
    report = None

    @classmethod
    def setUpClass(cls):
        if not fuzz_bodies():
            raise unittest.SkipTest('no cached bodies to fuzz')
        cls.workdir = tempfile.mkdtemp(prefix='fuzz-', dir=R.SCRATCH)
        cls.results, cls.report = run_child(cls.workdir)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(getattr(cls, 'workdir', ''), ignore_errors=True)
        r = cls.results or {}
        if r:
            fams = {}
            for x in r.values():
                fams.setdefault(x['family'], {}).setdefault(x['outcome'], 0)
                fams[x['family']][x['outcome']] += 1
            sys.stderr.write('\n[test_real_fuzz] %d cases over %d bodies in %.0f s: %s\n'
                             % (len(r), len({(x['cell'], x['i']) for x in r.values()}),
                                cls.report['wall_s'], json.dumps(fams, sort_keys=True)))

    def of(self, *families):
        return [r for r in self.results.values() if r['family'] in families]

    def test_no_hang_crash_or_runaway_allocation(self):
        rep = self.report
        self.assertIsNone(rep['killed'], 'the fuzz child was killed (%s) during %r'
                          % (rep['killed'], rep['in_progress']))
        self.assertEqual(0, rep['returncode'], 'the fuzz child crashed: %s'
                         % rep['stderr'])
        self.assertIsNotNone(rep['done'])
        self.assertEqual(rep['done']['done'], len(self.results))
        self.assertGreater(len(self.results), 1000)
        slow, fat, grow = [], [], []
        for r in self.results.values():
            if r['time_s'] > max(CASE_TIME[0], CASE_TIME[1] * r['base_time_s']):
                slow.append('%.3fs (base %.4fs) %s' % (r['time_s'], r['base_time_s'], _fmt(r)))
            if r['peak'] > max(CASE_PEAK[0], CASE_PEAK[1] * r['base_peak']):
                fat.append('%d B (base %d B) %s' % (r['peak'], r['base_peak'], _fmt(r)))
            if r['rss_growth'] > CASE_RSS_GROWTH:
                grow.append('%d MB %s' % (r['rss_growth'] >> 20, _fmt(r)))
        self.assertEqual([], slow[:10], '%d slow cases' % len(slow))
        self.assertEqual([], fat[:10], '%d cases over the allocation bound' % len(fat))
        self.assertEqual([], grow[:10], '%d cases grew the RSS' % len(grow))
        self.assertLess(rep['done']['max_rss'], RSS_CAP_BYTES)

    def test_corruption_is_a_format_error(self):
        """Truncation, a dropped or duplicated record, a count off by one, a huge or
        overflowing number, invalid UTF-8, CR LF: GraphletFormatError, every time."""
        families = [f for f in MUST if f not in PRODUCT_BUGS]
        bad = [_fmt(r) for r in self.of(*families) if r['outcome'] != 'format_error']
        self.assertEqual([], bad[:15], '%d of %d cases' % (len(bad), len(self.of(*families))))
        for f in families:
            self.assertTrue(self.of(f), 'no case of family %s' % f)

    def test_a_cut_document_is_never_a_graphlet(self):
        """Every proper prefix of every body fails: the Z record and the A counts make a
        truncated document detectable at any cut."""
        cuts = self.of('truncate_lines', 'truncate_bytes')
        self.assertGreater(len(cuts), 200)
        self.assertEqual([], [_fmt(r) for r in cuts if r['outcome'] == 'parsed'])

    def test_other_damage_parses_canonically_or_is_a_format_error(self):
        bad, parsed = [], 0
        for r in self.of(*MAY):
            if r['outcome'] == 'parsed':
                parsed += 1
                if r['roundtrip'] is not True:
                    bad.append('round trip %r: %s' % (r['roundtrip'], _fmt(r)))
            elif r['outcome'] != 'format_error':
                bad.append(_fmt(r))
        self.assertEqual([], bad[:15], '%d bad of %d' % (len(bad), len(self.of(*MAY))))

    def test_errors_name_a_line(self):
        """A format error names where: 'line N: ...' for a record-level defect; only a
        document-level defect (a cut document, the encoding) says no line."""
        bad = []
        for r in self.results.values():
            if r['outcome'] != 'format_error':
                continue
            if not re.match(r'line \d+: |truncated document|not UTF-8', r['msg']):
                bad.append(_fmt(r))
        self.assertEqual([], bad[:10], '%d without a line' % len(bad))

    def test_format_errors_are_short(self):
        """A message repeats at most 32 characters of a token: a 5,000-digit token must
        not make its own error huge (the tool layer would have to hide it)."""
        bad = [_fmt(r) for r in self.results.values()
               if r['outcome'] == 'format_error' and len(r['msg']) > MAX_MESSAGE]
        self.assertEqual([], bad[:5], '%d long messages' % len(bad))

    def test_run_errors_name_their_line(self):
        """The R-record range checks of parser._resolve_arm name the R record's line and
        the arm ('run 0' exists on both arms). Fixed product bug BUG-5: they said line
        0 and no arm. Repro of the old answer: swap two R records of
        sra_rand50_03__no_sequences (lines 11 and 16) -> 'run 0 ends at 25 outside its
        anchor 0 [0, 8]'."""
        runs = [r for r in self.results.values() if r['outcome'] == 'format_error'
                and RUN_RANGE_ERROR.match(r['msg'])]
        self.assertTrue(runs, 'no case reached the R range checks')
        bad = [_fmt(r) for r in runs
               if not re.match(r'line [1-9]\d*: (left|right) arm: run \d+ ends ', r['msg'])]
        self.assertEqual([], bad[:5], '%d of %d cases' % (len(bad), len(runs)))

    # the former product bugs: the cases are kept and named

    def test_huge_integer_tokens_are_a_format_error(self):
        # fixed BUG-1: parse('Z ' + '9' * 5000 ...) raised ValueError('Exceeds the limit
        # (4300 digits) ...'), not GraphletFormatError
        cases = self.of('huge_digits')
        self.assertTrue(cases)
        bad = [_fmt(r) for r in cases if r['outcome'] != 'format_error'
               or len(r['msg']) > MAX_MESSAGE]
        self.assertEqual([], bad[:5], '%d of %d' % (len(bad), len(cases)))

    def test_lone_surrogates_are_a_format_error(self):
        # fixed BUG-2: parse(<str with '\udc80' in an L name>) raised UnicodeEncodeError,
        # not GraphletFormatError; the error names the line of the surrogate
        cases = self.of('lone_surrogate')
        self.assertTrue(cases)
        bad = [_fmt(r) for r in cases if r['outcome'] != 'format_error'
               or not re.match(r'line [1-9]\d*: not UTF-8', r['msg'])]
        self.assertEqual([], bad[:5], '%d of %d' % (len(bad), len(cases)))

    def test_the_product_bug_cases_are_still_bounded_and_structured(self):
        """The former product-bug families fail fast, small, and as GraphletFormatError
        (which the tool layer turns into format_error)."""
        for r in self.of('huge_digits', 'lone_surrogate'):
            self.assertEqual('format_error', r['outcome'], _fmt(r))


# ------------------------------------------------------------------ the tool layer

class StubClient(TraverseClient):
    """The transport replaced by a canned (real, mutated) response: the tool layer's
    handling of the body is what is under test."""

    def __init__(self, response):
        super().__init__('stub', 0)
        self.response = response

    def traverse_raw(self, request):
        return copy.deepcopy(self.response)

    def capabilities(self):
        return dict(self.response.get('capabilities') or {})


class _Deadline:
    """SIGALRM-based guard for one in-process call (pure-Python loops only)."""

    def __init__(self, seconds):
        self.seconds = seconds

    def _fire(self, *a):
        raise TimeoutError('the call ran longer than %.0f s' % self.seconds)

    def __enter__(self):
        self.old = signal.signal(signal.SIGALRM, self._fire)
        signal.setitimer(signal.ITIMER_REAL, self.seconds)

    def __exit__(self, *exc):
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, self.old)
        return False


def tool_cases(per_family=3, families=FAMILIES, max_bodies=18):
    """A deterministic sample of (index, cell, i, family, desc, payload)."""
    bodies = fuzz_bodies()[:max_bodies] if max_bodies else fuzz_bodies()
    for index, cell, i, text in bodies:
        for family in families:
            muts = mutations(text, family, '%s/%s/%d' % (index, cell, i))
            for payload, desc in muts[:per_family]:
                yield index, cell, i, family, desc, payload


@R.skip_unless_cache()
class TestToolLayerFuzz(unittest.TestCase):
    """Corrupt real bodies through the agent's doors: traverse_fetch, graphlet_load and
    GraphletStore.put. A structured format_error (never a raise), nothing stored, no
    file content echoed, bounded time."""

    def setUp(self):
        self.spool = tempfile.mkdtemp(prefix='toolfuzz-', dir=R.SCRATCH)
        self.store = GraphletStore(self.spool)
        self.exports = os.path.join(self.spool, 'exports')
        os.makedirs(self.exports, exist_ok=True)

    def tearDown(self):
        shutil.rmtree(self.spool, ignore_errors=True)

    def state(self):
        return (sorted(e['handle'] for e in self.store.list()),
                sorted(os.listdir(os.path.join(self.spool, 'bodies'))))

    def fetch(self, index, cell, i, payload, *, fix_counts=True):
        """traverse_fetch returning the cell's real response with results[i]'s body
        replaced (graphlet_bytes / graphlet_lines recomputed when fix_counts, so the
        parse itself is reached)."""
        c = R.load_cell(index, cell)
        resp = copy.deepcopy(c.graphlet_response)
        r = resp['results'][i]
        resp['results'] = [r]
        text = payload.decode('utf-8', 'surrogateescape') if isinstance(payload, bytes) \
            else payload
        r['graphlet'] = text
        if fix_counts:
            # a lone surrogate has no UTF-8: count it as the server would have sent it
            r['graphlet_bytes'] = len(text.encode('utf-8', 'surrogatepass'))
            r['graphlet_lines'] = text.count('\n')
        tools = GraphletTools(self.store, {'stub': StubClient(resp)})
        before = self.state()
        with _Deadline(30):
            t0 = time.perf_counter()
            out = tools.traverse_fetch('stub', seed={'sequence': c.seed_sequence(i)})
            dt = time.perf_counter() - t0
        self.assertLess(dt, 10.0)
        self.assertLessEqual(len(json.dumps(out, separators=(',', ':'), ensure_ascii=False)
                                 .encode('utf-8', 'surrogatepass')), 2048)
        if 'error' in out:
            self.assertEqual(before, self.state(), 'a failed fetch stored something')
        return out

    def check_family(self, families, want=('format_error',), allow_parse=False):
        n, bad = 0, []
        for index, cell, i, family, desc, payload in tool_cases(families=families):
            out = self.fetch(index, cell, i, payload)
            n += 1
            if 'error' in out:
                if out['error'] not in want:
                    bad.append('%s/%d %s [%s]: %r' % (cell, i, family, desc, out))
            elif not allow_parse:
                bad.append('%s/%d %s [%s]: stored as %s' % (cell, i, family, desc,
                                                           out.get('handle')))
            else:
                self.store.free(out['handle'])
        self.assertGreater(n, 0)
        self.assertEqual([], bad[:3], '%d of %d' % (len(bad), n))

    def test_fetch_of_a_corrupt_body_is_a_format_error(self):
        self.check_family([f for f in MUST if f not in PRODUCT_BUGS and f != 'bad_utf8'])

    def test_fetch_of_damaged_records_parses_or_is_a_format_error(self):
        self.check_family(MAY, allow_parse=True)

    def test_transport_truncation_is_detected_by_the_envelope_counts(self):
        """A body cut in transport, its graphlet_bytes / graphlet_lines as sent: refused
        before anything is stored."""
        n = 0
        for index, cell, i, text in fuzz_bodies()[:20]:
            for payload, desc in mutations(text, 'truncate_lines', cell)[-3:]:
                out = self.fetch(index, cell, i, payload, fix_counts=False)
                self.assertEqual('format_error', out.get('error'), (cell, desc, out))
                n += 1
        self.assertGreater(n, 10)

    def test_store_put_rejects_and_stores_nothing(self):
        n = 0
        for index, cell, i, family, desc, payload in tool_cases(
                families=('truncate_lines', 'drop_line', 'count_flip', 'huge_range', 'crlf'),
                per_family=2):
            c = R.load_cell(index, cell)
            resp = copy.deepcopy(c.graphlet_response)
            r = resp['results'][i]
            r['graphlet'] = payload
            r['graphlet_bytes'] = len(payload.encode('utf-8'))
            r['graphlet_lines'] = payload.count('\n')
            before = self.state()
            with self.assertRaises(GraphletFormatError, msg=(cell, family, desc)):
                self.store.put(resp, c.request, i)
            self.assertEqual(before, self.state())
            n += 1
        self.assertGreater(n, 20)

    def test_load_of_a_corrupt_file_is_a_format_error(self):
        """Saved .mgt files (H, J, body) damaged on disk: format_error naming a line and
        nothing of the file's content (it may be anything the agent named)."""
        tools = GraphletTools(self.store, {})
        n, bad = 0, []
        for index, cell, i, text in fuzz_bodies()[:18]:
            g = R.load_cell(index, cell).graphlet(i)
            saved = dump(g, envelope=True)
            names = [l.name for l in g.labels if len(l.name) >= 6]
            variants = []
            for family in ('truncate_lines', 'truncate_bytes', 'drop_line', 'dup_line',
                           'count_flip', 'huge_range', 'bad_utf8', 'crlf'):
                variants += [(family, d, p) for p, d in mutations(saved, family, cell)[:2]]
            lines = _lines(saved)
            for j in ('J {', 'J []', 'J null', 'J "x"', 'J {"view": 7}'):
                variants.append(('j_corrupt', j, _join(lines[:1] + [j] + lines[2:])))
            for family, desc, payload in variants:
                name = 'f%d.mgt' % n
                with open(os.path.join(self.exports, name), 'wb') as f:
                    f.write(payload if isinstance(payload, bytes) else payload.encode())
                before = self.state()
                with _Deadline(30):
                    out = tools.graphlet_load(name)
                n += 1
                if out.get('error') != 'format_error':
                    if family == 'j_corrupt' and 'handle' in out:
                        tools.graphlet_free(out['handle'])
                        continue
                    bad.append('%s %s [%s]: %r' % (cell, family, desc, out))
                    continue
                self.assertEqual(before, self.state())
                self.assertRegex(out['message'], r'line \d+|not UTF-8')
                for nm in names:
                    self.assertNotIn(nm, out['message'])
        self.assertGreater(n, 100)
        self.assertEqual([], bad[:10], '%d bad' % len(bad))

    def test_invalid_utf8_through_the_tools(self):
        """Bytes that are not UTF-8: a file -> format_error; a transport that decodes
        them with surrogateescape -> a structured error, nothing stored."""
        tools = GraphletTools(self.store, {})
        n = 0
        for index, cell, i, family, desc, payload in tool_cases(families=('bad_utf8',),
                                                                per_family=5):
            with open(os.path.join(self.exports, 'u.mgt'), 'wb') as f:
                f.write(payload)
            out = tools.graphlet_load('u.mgt')
            self.assertEqual('format_error', out['error'], (cell, desc, out))
            out = self.fetch(index, cell, i, payload)
            self.assertIn('error', out, (cell, desc))
            n += 1
        self.assertGreater(n, 5)

    def test_fetch_reports_huge_integers_as_format_error(self):
        # fixed BUG-1 through the tool layer: traverse_fetch answered {'error':
        # 'bad_argument', 'message': 'Exceeds the limit (4300 digits) ...'} -- it blamed
        # the agent's arguments for a corrupt server body
        self.check_family(['huge_digits'])

    def test_fetch_reports_lone_surrogates_as_format_error(self):
        # fixed BUG-2 through the tool layer: traverse_fetch answered bad_argument
        # ("'utf-8' codec can't encode character '\udc80' ...") for a body whose JSON
        # string carries a lone surrogate escape
        self.check_family(['lone_surrogate'])

    def test_the_product_bug_bodies_are_still_structured_and_not_stored(self):
        """The former BUG-1/BUG-2 bodies: format_error (never bad_argument, never a
        raise, never result_too_large), nothing stored."""
        n = 0
        for index, cell, i, family, desc, payload in tool_cases(
                families=('huge_digits', 'lone_surrogate'), per_family=2):
            out = self.fetch(index, cell, i, payload)
            self.assertEqual('format_error', out.get('error'), (cell, desc, out))
            n += 1
        self.assertGreater(n, 5)

    def test_a_huge_token_keeps_its_error_class(self):
        """FORMAT-ERROR-ECHO: a corrupt body whose bad token is thousands of characters
        long is still a format_error under the 2 KB ceiling (the message echoes 32
        characters of it), not a result_too_large that hides the error class."""
        n = 0
        for index, cell, i, text in fuzz_bodies()[:6]:
            for junk in ('x' * 3000, '9' * 5000):
                body = re.sub(r'^Z \d+$', 'Z ' + junk, text, flags=re.M)
                self.assertNotEqual(text, body)
                out = self.fetch(index, cell, i, body)
                self.assertEqual('format_error', out.get('error'), (cell, out))
                self.assertLess(len(out['message']), MAX_MESSAGE, out)
                n += 1
        self.assertGreater(n, 0)


if __name__ == '__main__':
    if len(sys.argv) == 3 and sys.argv[1] == '--fuzz-child':
        child_main(sys.argv[2])
    else:
        unittest.main()
