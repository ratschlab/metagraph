"""The rules of GET /traverse/capabilities and GET /capabilities are references to the sections of
docs/SPEC-labeled-traversal-core.md that state them (its §10.3, "The rules are references"),
checked on the capabilities documents of the mini index's fixture servers
(data/traverse/pattern/, written by scripts/traversal/pattern_fixtures.py) against the SPEC's
text, never against the server's sources:

  - the table of §10.3 names every rule field, its route and its section, and every fixture
    document states exactly that reference in that field (the probe's fields on
    /traverse/capabilities, the others on both routes);
  - every referenced section exists and holds each named part of the reference;
  - no other long text is left in the documents (a rule written out again would be);
  - each referenced section states the rule's substance: the phrases below, which the rules'
    texts carried and the tests checked there before (normalised: no Markdown emphasis or code
    marks, ASCII dashes, single spaces, lower case).
"""

import json
import os
import re
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data', 'traverse', 'pattern')
REPO = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
SPEC = os.path.join(REPO, 'docs', 'SPEC-labeled-traversal-core.md')
SPEC_NAME = 'SPEC-labeled-traversal-core.md'

# the longest string a capabilities document may hold: a reference, an id pattern, a route, a
# path of a fixture server
LONGEST_TEXT = 100

# What each referenced section states, phrase by phrase (normalised as normalise() does)
SUBSTANCE = {
    'work_bound': [
        'deterministic logical work, not measured decode effort',
        'successor enumerations (4)', '8 per key, 1 per entry', 'dependency rows',
        'whatever the decode shared or cached',
        'the most its seed charged between two comparisons',
        'at least w units over budget', 'the deadline is read before every head',
        'physical decode counters',
    ],
    'decode_cache.rule': [
        '--traverse-path-cache-mb', 'in two generations', 'never depends on the cache',
        'every row it was asked for', 'the first 8 rows after each of them on its path',
        'multiple of 16', '4,096 bytes', 'kept flat', 'min(budget / 4, 64 mib)',
        "without a memory budget it is the request's, kept from seed to seed",
    ],
    'coordinates.rule': [
        'coordinates_reason', '0-based, half-open, forward strand',
        '[c + k - l, c + k) on the right arm, [c, c + l) on the left',
        "a record's last k - 1 bases share their numbers with the next record's first k - 1 "
        'positions',
        'occurrences_total', 'chains_ended', 'lower_bound: true', 'one strand per record',
        'unique headers',
    ],
    'coordinates.output_bound': [
        'at most min(chains, max_coordinate_occurrences) per run and seed label',
        'which a request without one runs under when it is not 0',
        'the coordinates read are within it', '0, their default',
        'a client sets a budget, keeps the cap or drops the coordinates',
    ],
    'deadline_check.rule': [
        'decoded in chunks, the deadline checked before each',
        'less than 1/64 of the time left', 'first chunk is at most 8 rows',
        'at most 4 times the previous one', 'less than a quarter of the time left',
        "a seed's own time budget never stops its validation", 'it stays null',
        'returns exactly what one read returns', 'every 16 graph steps', 'every 64 kib',
        'not chunked', 'decodes every read in one piece',
        'so up to poll_stride - 1 heads after it', 'the poll before each seed',
    ],
    'resolve.time_budget.rule': [
        'opt-in', 'above the finalisation reserve', 'limits.clamped',
        'the work stops at the budget less the reserve',
        'a stop answers the resolve of a query prefix', 'remainder_from_bp', 'selection: null',
        'the answer is built within the reserve', 'not polled',
        "the mapping of the query's k-mers", 'select_seeds',
    ],
    'attempts.bound': [
        "which starts when the server read the request's header",
        'bound_ms - max(allowance_ms / 2, reserve_ms)',
        'its first poll that reads the clock after it', 'one in poll_stride',
        'the bound itself is compared only at those delivery checks',
        'past the bound the attempt runs on until its next delivery check',
        'a run whose length has no stated bound', '503 with usage',
        'a lower walk-until in force only between two such polls',
    ],
    'attempts.not_after': [
        'cannot start subsequently', 'not_after_ms + attempts.clock_skew_allowance_ms',
        'the check is strict', 'the refusals are made in one order, under one lock',
        'instance_mismatch', 'the expiry last',
        'so an expired 409 always means that no attempt with that id existed on that '
        'server_instance when the request was judged',
        'apart from what it runs past its bound up to its next delivery check',
        'a later get /traverse/attempt/{id} answers 404',
    ],
    'attempts.instance': [
        'compared as given', 'checked first, ahead of a duplicate id and of not_after_ms',
        'a delayed copy of a cancelled or of a finished request would otherwise run there',
        'without attempt_id it is a 400 naming the field',
        'nor does a server below feature level 5', 'exit status 1',
    ],
    'attempts.suppression': [
        'tombstoned', 'held at least retention_s from the cancel',
        "until this server's clock reads not_after_ms + clock_skew_allowance_ms",
        'at most tombstone_max_s from the cancel', 'never shorten it',
        'while either clock says it is live', 'while it reads at most suppressed_until_ms',
        'covers_admission (true iff a not_after_ms was judged against',
        'covers_admission_reason (no_not_after_ms | beyond_tombstone_max)',
        "for the refused request's 409 its own alone, absent when it has none",
        'tombstones_full', 'no_suppression', 'what a clock step can do',
    ],
    'attempts.release_rule': [
        'covers_admission: true from the same server_instance', 'exactly that not_after_ms',
        'never released early on a tombstone', 'the 409 refusing a copy of the request whose',
        'not_after_ms + clock_skew_allowance_ms + bound_ms',
        "a finished attempt's id stays refused (409", 'they are never dropped early',
        "the hold is the process's, in memory", 'replay-safe',
        'pins the instance, as for a tombstone',
        'assumes that no copy of the request arrives after the attempt left retention',
        'with retention_s 0 nothing is kept',
        'a pinned copy is still refused by a restarted process',
        'it settles nothing', 'server_time_ms - not_after_ms is the step it survives',
        'the clock release takes the attempt as stopped',
        'apart from what it runs past it until its next delivery check',
        'a run whose length has no stated bound',
        'keeps polling while it answers so, or adds a margin of its own',
        'a 409 instance_mismatch release nothing', 'another server that serves the same ledger',
        "that the attempt's run past its bound has ended by then",
        'at most tombstone_max_s after it finished or after its latest refused copy',
        'the finishes and refused copies within tombstone_max_s',
    ],
    'attempts.delivery_reserve.rule': [
        'reserve_ms = 1.25 x ((t + e) / (compress_mbps x 1000) + e / (build_mbps x 1000)) '
        '+ stop_ms',
        'coordinate_account_per_text_byte', "published at every level's end",
        'the slowest rate', "and the attempt's own seeds' where more conservative",
        'on a rest of at least 1 mib', 'measured_stop_ms', 'when a 503 becomes likely',
        'writes more text than its account / account_per_text_byte',
    ],
    'attempts.delivery_reserve.calibration': [
        'starting estimates', '33.5', '1,344', '58.4', '352 ms', '1,001 ms', '1,699 ms',
        '460-670 mb/s', '13.6-51 mb/s', '115-129', '4 x build_mbps',
    ],
}


def normalise(text):
    """Markdown to the words it says: no emphasis or code marks, ASCII dashes and times, single
    spaces, lower case."""
    text = text.replace('`', '').replace('*', '')
    for a, b in (('−', '-'), ('–', '-'), ('—', '-'), ('×', 'x'),
                 (' ', ' '), (' ', ' ')):
        text = text.replace(a, b)
    return re.sub(r'\s+', ' ', text).strip().lower()


def sections(spec):
    """{'N' or 'N.M': its text}, a section of level 2 with its subsections."""
    out = {}
    heads = [(m.start(), m.group(1), m.group(2))
             for m in re.finditer(r'^(#{2,3}) (\d+(?:\.\d+)?)\.? ', spec, re.M)]
    for i, (start, level, number) in enumerate(heads):
        end = len(spec)
        for start2, level2, _ in heads[i + 1:]:
            if len(level2) <= len(level):
                end = start2
                break
        out[number] = spec[start:end]
    return out


def reference_table(spec):
    """[(field, route, section, name)] of §10.3's table "The rules are references"."""
    m = re.search(r'\*\*The rules are references\.\*\*.*?\n\n(.*?)\n\n', spec, re.S)
    assert m, 'the table of §10.3, "The rules are references"'
    rows = []
    for line in m.group(1).strip().split('\n')[2:]:
        field, route, where = [c.strip() for c in line.strip().strip('|').split('|')]
        number, _, name = where.partition(',')
        rows.append((field.strip('`'), route, number.strip().lstrip('§'),
                     name.strip().replace('`', '') or None))
    return rows


def reference(number, name):
    return f'{SPEC_NAME} section {number}' + (f', {name}' if name else '')


def at(doc, field):
    for key in field.split('.'):
        if not isinstance(doc, dict) or key not in doc:
            return None
        doc = doc[key]
    return doc


def strings(value):
    if isinstance(value, dict):
        for v in value.values():
            yield from strings(v)
    elif isinstance(value, list):
        for v in value:
            yield from strings(v)
    elif isinstance(value, str):
        yield value


class TestTraverseCapabilitiesReferences(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        with open(os.path.join(DATA, 'index.json'), encoding='utf-8') as f:
            fixtures = json.load(f)['fixtures']
        cls.documents = {}
        for name, f in fixtures.items():
            route = f['path'].split('?')[0]
            if f['method'] == 'GET' and f['status'] == 200 \
                    and route in ('/traverse/capabilities', '/capabilities'):
                with open(os.path.join(DATA, name, 'answer.json'), encoding='utf-8') as fh:
                    cls.documents[name] = (route, json.load(fh))
        cls.spec = None
        if os.path.isfile(SPEC):
            with open(SPEC, encoding='utf-8') as f:
                cls.spec = f.read()

    def setUp(self):
        if self.spec is None:
            self.skipTest('docs/SPEC-labeled-traversal-core.md is not in this checkout')

    def check_documents(self, documents, table):
        """The fixture documents state the table's references, and no other long text."""
        for name, (route, doc) in documents.items():
            for field, where, number, part in table:
                value = at(doc, field)
                if where == 'probe' and route != '/traverse/capabilities':
                    self.assertIsNone(value, f'{name}: {field} only on the probe')
                    continue
                self.assertEqual(reference(number, part), value, f'{name}: {field}')
            for text in strings(doc):
                self.assertLessEqual(len(text), LONGEST_TEXT, f'{name}: {text[:80]!r}...')

    def check_sections(self, table, by_number):
        """Every referenced section exists, holds each named part and states the substance."""
        for field, _, number, part in table:
            self.assertIn(number, by_number, f'{field}: section {number}')
            text = normalise(by_number[number])
            for piece in (part or '').split(': '):
                if piece:
                    self.assertIn(normalise(piece), text, f'{field}: section {number}, {piece}')
            for phrase in SUBSTANCE[field]:
                self.assertIn(phrase, text, f'{field}: section {number}')

    def test_the_documents_name_the_sections(self):
        table = reference_table(self.spec)
        self.assertEqual(set(SUBSTANCE), {field for field, _, _, _ in table})
        self.assertEqual({'both', 'probe'}, {route for _, route, _, _ in table})
        # every capabilities document of every fixture server, both routes
        self.assertEqual({'/traverse/capabilities', '/capabilities'},
                         {route for route, _ in self.documents.values()})
        self.check_documents(self.documents, table)
        self.check_sections(table, sections(self.spec))

    def test_a_reference_that_names_nothing_fails(self):
        """The checks are not vacuous: a reference to a section that does not exist, a part the
        section does not hold, a rule the section does not state, a field that is not the
        reference, and a text written out again each fail."""
        table = reference_table(self.spec)
        by_number = sections(self.spec)
        broken = [
            [(f, r, '6.99' if f == 'work_bound' else n, p) for f, r, n, p in table],
            [(f, r, n, 'chunked windows' if f == 'deadline_check.rule' else p)
             for f, r, n, p in table],
            [(f, r, '4.5' if f == 'attempts.release_rule' else n,
              None if f == 'attempts.release_rule' else p) for f, r, n, p in table],
        ]
        for rows in broken:
            with self.assertRaises(AssertionError):
                self.check_sections(rows, by_number)
        name, (route, doc) = next((n, d) for n, d in self.documents.items()
                                  if d[0] == '/traverse/capabilities')
        changed = json.loads(json.dumps(doc))
        changed['attempts']['bound'] = reference('6.8', 'the delivery reserve')
        with self.assertRaises(AssertionError):
            self.check_documents({name: (route, changed)}, table)
        changed = json.loads(json.dumps(doc))
        changed['work_bound'] += '; ' + 'x' * LONGEST_TEXT
        with self.assertRaises(AssertionError):
            self.check_documents({name: (route, changed)}, table)
        changed = json.loads(json.dumps(doc))
        changed['attempts']['id_pattern'] = 'y' * (LONGEST_TEXT + 1)
        with self.assertRaises(AssertionError):
            self.check_documents({name: (route, changed)}, table)


if __name__ == '__main__':
    unittest.main()
