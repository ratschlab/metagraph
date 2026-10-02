"""MGT v1 reading and writing: the whole-document fixtures (§2.7), the canonical form
(§2.3), truncation and malformed input, the envelope (J, from_response, save/load)."""

import json
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import traverse_testlib as T  # noqa: E402

from metagraph.traverse import (  # noqa: E402
    GraphletFormatError, MissingEnvelope, dump, from_response, is_canonical, load, parse,
)
from metagraph.traverse.export import normalize_result  # noqa: E402


class TestDocuments(unittest.TestCase):
    """For every fixture: dump(parse(x)) == x byte-exactly, the detail: graphlet
    response carries that body, and to_json() reproduces the detail: full result."""

    def test_all_fixtures_are_present(self):
        names = sorted(f[:-4] for f in os.listdir(T.DOCS) if f.endswith('.mgt'))
        self.assertEqual(sorted(T.DOCUMENTS), names)

    def test_dump_is_identity_and_canonical(self):
        for name in T.DOCUMENTS:
            with self.subTest(name):
                text = T.doc_text(name)
                g = parse(text)
                self.assertEqual(text, g.dump())
                self.assertTrue(is_canonical(text))
                self.assertEqual(text, dump(parse(text.encode('utf-8'))))

    def test_response_carries_the_body(self):
        for name in T.DOCUMENTS:
            with self.subTest(name):
                text = T.doc_text(name)
                resp = T.doc_json(name, 'graphlet')
                r = resp['results'][0]
                self.assertEqual(text, r['graphlet'])
                self.assertEqual(len(text.encode('utf-8')), r['graphlet_bytes'])
                self.assertEqual(text.count('\n'), r['graphlet_lines'])
                self.assertEqual('Z %d' % r['graphlet_lines'], text.split('\n')[-2])
                self.assertEqual('graphlet', resp['strategy']['output']['detail'])

    def test_to_json_reproduces_detail_full(self):
        for name in T.DOCUMENTS:
            with self.subTest(name):
                full = T.doc_json(name, 'full')['results'][0]
                mine = normalize_result(T.graphlet(name).to_json())
                theirs = normalize_result(full)
                self.assertEqual(theirs, mine, '\n'.join(T.json_diff(mine, theirs)[:20]))

    def test_check_rules_finds_no_divergence(self):
        for name in T.DOCUMENTS:
            with self.subTest(name):
                full = T.doc_json(name, 'full')['results'][0]
                self.assertEqual([], T.graphlet(name).check_rules(full))

    def test_a_counts_match_the_body(self):
        for name in T.DOCUMENTS:
            g = T.body(name)
            for arm in g.arms.values():
                self.assertEqual(len(arm.segments), arm.counts.segments)
                self.assertEqual(len(arm.runs), arm.counts.runs)
                self.assertEqual(sum(s.length_bp for s in arm.segments), arm.counts.bases)

    def test_label_end_events_never_outnumber_runs(self):
        for name in T.DOCUMENTS:
            full = T.doc_json(name, 'full')['results'][0]
            for arm in full['arms'].values():
                ends = sum(1 for s in arm['segments'] for e in s['events']
                           if e['type'] == 'label_end')
                self.assertLessEqual(ends, len(arm['runs']))


class TestEnvelope(unittest.TestCase):
    def test_body_only_raises_missing_envelope(self):
        g = T.body('fork')
        self.assertFalse(g.has_envelope)
        with self.assertRaises(MissingEnvelope):
            g.to_json()
        with self.assertRaises(MissingEnvelope):
            g.summary()
        with self.assertRaises(MissingEnvelope):
            g.next_request('right', [0])
        # the body-level queries still work
        self.assertEqual(2, len(g.walks('right')))

    def test_save_load_round_trip_with_j(self):
        g = T.graphlet('merge')
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, 'merge.mgt')
            g.save(path)
            text = T.read(path)
            lines = text.split('\n')
            self.assertTrue(lines[1].startswith('J {'))
            self.assertEqual('Z %d' % (len(lines) - 1), lines[-2])
            j = json.loads(lines[1][2:])
            self.assertEqual(1, len(j['results']))
            self.assertNotIn('graphlet', j['results'][0])
            self.assertEqual(j, json.loads(json.dumps(j, sort_keys=True)))
            back = load(path)
            self.assertTrue(back.has_envelope)
            self.assertEqual(text, back.dump())               # identity with J
            self.assertEqual(T.doc_text('merge'), back.dump(envelope=False))
            self.assertEqual(normalize_result(g.to_json()), normalize_result(back.to_json()))
            self.assertTrue(is_canonical(text))

    def test_j_must_follow_h_once(self):
        g = T.graphlet('linear')
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, 'x.mgt')
            g.save(path)
            lines = T.read(path).split('\n')
        twice = lines[:2] + [lines[1]] + lines[2:-2] + ['Z %d' % (len(lines))] + ['']
        with self.assertRaises(GraphletFormatError):
            parse('\n'.join(twice))
        late = [lines[0], lines[2], lines[1]] + lines[3:]
        with self.assertRaises(GraphletFormatError):
            parse('\n'.join(late))

    def test_from_response_rejects_a_transport_cut(self):
        resp = T.doc_json('fork', 'graphlet')
        r = dict(resp['results'][0])
        r['graphlet_bytes'] += 1
        with self.assertRaises(GraphletFormatError):
            from_response(r, resp)
        with self.assertRaises(ValueError):
            from_response({'seed': {}, 'error': 'no label carries every k-mer'}, resp)


class TestTruncation(unittest.TestCase):
    """A truncated body raises: never a shallower graphlet (§14)."""

    def test_every_prefix_raises(self):
        for name in ('merge', 'annotate', 'limits'):
            lines = T.doc_text(name).split('\n')[:-1]
            for k in range(1, len(lines)):
                with self.subTest(name=name, lines=k):
                    with self.assertRaises(GraphletFormatError):
                        parse('\n'.join(lines[:k]) + '\n')

    def test_missing_final_line_feed(self):
        with self.assertRaises(GraphletFormatError):
            parse(T.doc_text('linear')[:-1])

    def test_dropping_one_record_raises(self):
        """Even with Z rewritten to match: the A counts, the bins (steps and label ends),
        the branch-event total and the continuation rule notice a lost record. A K
        record is the exception -- it is optional, and check_rules() reports the outcome
        that its loss leaves inconsistent."""
        # not counted anywhere: the label dictionary, the resource stop, revisit and
        # switch events (a switch into a split child repeats what the child's entry set
        # and its run's from_label already state), a V record beyond an evidence boundary
        # (the cap that bounds them is in the envelope), and K records -- which
        # check_rules() follows
        def uncounted(line, evidence_cut):
            f = line.split(' ')
            return (f[0] in ('L', 'Q', 'K') or (f[0] == 'E' and f[2] in ('s', 'v'))
                    or (f[0] == 'V' and evidence_cut))

        for name in ('merge', 'limits', 'annotate', 'switch_chain', 'reentry', 'caps'):
            lines = T.doc_text(name).split('\n')[:-1]
            evidence_cut = False
            for i in range(1, len(lines) - 1):
                if lines[i].startswith('A '):
                    evidence_cut = lines[i].split(' ')[13] != '*'
                cut = lines[:i] + lines[i + 1:]
                cut[-1] = 'Z %d' % len(cut)
                text = '\n'.join(cut) + '\n'
                with self.subTest(name=name, record=lines[i][:30]):
                    if uncounted(lines[i], evidence_cut):
                        try:                 # may or may not be noticed
                            parse(text)
                        except GraphletFormatError:
                            pass
                        continue
                    with self.assertRaises(GraphletFormatError):
                        parse(text)

    def test_check_rules_notices_a_lost_limitation(self):
        for name, prefix in (('merge', 'K r scope'), ('limits', 'K r branch_events'),
                             ('limits', 'K * trace'), ('annotate', 'K r label_lists'),
                             ('caps', 'K r walk_domain')):
            lines = T.doc_text(name).split('\n')[:-1]
            i = next(j for j, l in enumerate(lines) if l.startswith(prefix))
            if name == 'annotate':
                # every label-class record of both arms keeps the lower bound: drop them all
                lines = [l for l in lines if not l.startswith(('K l label_lists', 'K l inexact',
                                                               'K r label_lists', 'K r inexact'))]
            else:
                lines = lines[:i] + lines[i + 1:]
            lines[-1] = 'Z %d' % len(lines)
            g = parse('\n'.join(lines) + '\n')
            with self.subTest(name=name, record=prefix):
                self.assertTrue(any(f['field'].startswith('outcome') for f in g.check_rules()))

    def test_a_counts_catch_a_corrupt_arm(self):
        t = T.doc_text('fork')
        with self.assertRaises(GraphletFormatError) as e:
            parse(_edit(t, ' * 0 * 3 3 2 1 0 20\n', ' * 0 * 3 3 2 1 0 21\n'))
        self.assertIn('A counts', str(e.exception))


def _edit(text, old, new, count=1):
    assert old in text, old
    return text.replace(old, new, count)


class TestFuzz(unittest.TestCase):
    def test_random_edits_parse_or_raise_a_format_error(self):
        """Any byte edited, deleted or inserted: a GraphletFormatError, or a document
        whose canonical dump is a fixed point -- never another exception."""
        import random
        rng = random.Random(5)
        alphabet = '0123456789 *.!+-,|:;=%ACGTXDRLmwyhcabsvu\n'
        for name in T.DOCUMENTS:
            t = T.doc_text(name)
            for _ in range(300):
                i = rng.randrange(len(t))
                op = rng.random()
                if op < 0.4:
                    m = t[:i] + rng.choice(alphabet) + t[i + 1:]
                elif op < 0.7:
                    m = t[:i] + t[i + 1:]
                else:
                    m = t[:i] + rng.choice(alphabet) + t[i:]
                try:
                    d = parse(m).dump()
                except GraphletFormatError:
                    continue
                self.assertEqual(d, parse(d).dump(), (name, i))


class TestCanonicalForm(unittest.TestCase):
    """A valid non-canonical document parses and dumps to its canonical form (§2.3)."""

    def assert_noncanonical(self, canonical, variant):
        self.assertNotEqual(canonical, variant)
        self.assertFalse(is_canonical(variant))
        self.assertEqual(canonical, parse(variant).dump())

    def test_explicit_values_where_star_does(self):
        t = T.doc_text('linear')
        # the right root's entry, entry_total and end written out (acc2 ends at 0)
        root = 'G * 0 40 * * * * * * TTACACCACGAGATCGCCTGGGGCTCTGACAGTTAGCATA'
        self.assert_noncanonical(t, _edit(t, root, root.replace(' * * * * * * ',
                                                                ' 0-2 3 0,2 * * * ')))

    def test_longer_setexpr_spelling(self):
        t = T.doc_text('merge')
        self.assert_noncanonical(t, _edit(t, 'G 0 20 16 !1 ', 'G 0 20 16 0,2-3 '))
        self.assert_noncanonical(t, _edit(t, 'G 0 20 16 1,3 ', 'G 0 20 16 !0,2 '))

    def test_maps_in_canonical_order(self):
        t = T.doc_text('switch_chain')
        self.assert_noncanonical(t, _edit(t, ' L:1,X:2', ' X:2,L:1'))
        t = T.doc_text('fork')
        self.assert_noncanonical(t, _edit(t, 'T X .', 'T X 1:0:0:0'))

    def test_shorter_front_coding_prefix(self):
        t = T.doc_text('linear')
        self.assert_noncanonical(t, _edit(t, 'L h 0 1 3 2', 'L h 0 1 0 acc2'))

    def test_first_base_only_on_split_children(self):
        t = T.doc_text('limits')
        self.assert_noncanonical(t, _edit(t, 'G * 0 3 * * * * 0 * *', 'G * 0 3 * * * * 0 G *'))
        with self.assertRaises(GraphletFormatError):          # a split child needs it
            parse(_edit(t, 'G 0 3 4 !2 * * * * C *', 'G 0 3 4 !2 * * * * * *'))

    def test_explicit_continuation_beside_bases(self):
        t = T.doc_text('fork')
        # the right arm's first leaf (path 0): its continuation crosses into the seed
        derived = parse(t).continuation('right', 0).sequence
        self.assertEqual(15, len(derived))
        right = t.index('\nA r ')
        c = t.index('\nC 15 0 0 1\n', right) + 1
        line = 'C 15 0 0 1\n'
        # equal to the derived spelling: valid, not canonical
        self.assert_noncanonical(t, t[:c] + 'C 15 0 0 1 %s\n' % derived + t[c + len(line):])
        # different: the writer's fallback, kept and preferred over the derivation
        g = parse(t[:c] + 'C 15 0 0 1 %s\n' % ('G' * 15) + t[c + len(line):])
        self.assertEqual('G' * 15, g.continuation('right', 0).sequence)
        self.assertTrue(is_canonical(g.dump()))

    def test_annotate_partition_written_out(self):
        t = T.doc_text('annotate')
        self.assert_noncanonical(t, _edit(t, 'G 1,2 36 10 ! 4 * * 0', 'G 1,2 36 10 ! 4 * .|. 0'))


class TestMalformed(unittest.TestCase):
    def assert_bad(self, text, msg=None):
        with self.assertRaises(GraphletFormatError) as e:
            parse(text)
        if msg:
            self.assertIn(msg, str(e.exception))
        return e.exception

    def test_version_and_tags(self):
        t = T.doc_text('linear')
        self.assert_bad(_edit(t, 'H mgt 1 ', 'H mgt 2 '), 'unsupported MGT version')
        self.assert_bad(_edit(t, ' walk ', ' natural '), 'orientation')
        self.assert_bad(_edit(t, 'O c c c i', 'O c c c z'), 'delivery')
        self.assert_bad(_edit(t, 'T * .', 'W * .'))

    def test_line_numbers(self):
        t = T.doc_text('linear')
        e = self.assert_bad(_edit(t, 'R 0 1 0 40 D 0 * * 0 0 0', 'R 0 1 0 40 D 0 * * 0 0 x'))
        self.assertEqual(18, e.line_no)

    def test_records_out_of_place(self):
        t = T.doc_text('linear')
        root = 'G * 0 40 * * * * * * GCACTCCCTATGAAAACGTCTGTCAGCTACGAGGAGACTA\n'
        self.assert_bad(_edit(t, root, root + 'P 0 1 1 0\n'), 'annotate mode only')
        t = T.doc_text('fork')
        self.assert_bad(_edit(t, 'T X .\nC 15 0 0 !\n', 'C 15 0 0 !\n'), 'without its T')

    def test_run_fields(self):
        t = T.doc_text('limits')
        self.assert_bad(_edit(t, 'R 1 1 0 5 B 0 * * 2 0 0 2', 'R 1 1 0 5 B 0 * * 2 0 0'),
                        'needed_budget')
        t = T.doc_text('switch_chain')
        self.assert_bad(_edit(t, 'R 0 0 0 30 Lw 0 * * * 0 0', 'R 0 0 0 30 Lw 0 * * 2 0 0'),
                        'structural_successors')
        # run 2 cannot follow run 3
        self.assert_bad(_edit(t, 'R 1 3 60 80 X 0 1:1 1 0 1 2', 'R 1 3 60 80 X 0 1:1 3 0 1 2'),
                        'prev_run')
        self.assert_bad(_edit(t, 'R 1 3 60 80 X 0 1:1 1 0 1 2', 'R 1 7 60 80 X 0 1:1 1 0 1 2'),
                        'label id')

    def test_structure(self):
        t = T.doc_text('merge')
        self.assert_bad(_edit(t, '0,2-3|1 1', '0,2-3 1'), 'partition')
        self.assert_bad(_edit(t, 'G * 0 20 * * * * 1 * ', 'G * 0 20 * * * * * * '), 'G.split')
        self.assert_bad(_edit(t, 'G 0 20 16 !1 ', 'G 0 20 16 !+1 '), 'does not fit its base')
        self.assert_bad(_edit(t, 'G 0 20 16 !1 ', 'G 0 20 16 !9 '), 'beyond the 4 L records')
        t = T.doc_text('limits')
        self.assert_bad(_edit(t, 'G 0 3 4 2 * * * * G *', 'G 0 3 4 2 * * * * * GACT'),
                        'bases on some segments only')

    def test_ranges_are_bounded_before_they_are_expanded(self):
        # a corrupt RANGES token of a few bytes is rejected against the label count before
        # it is expanded: no allocation in proportion to the range (§14, a local failure)
        import tracemalloc
        t = T.doc_text('fork')
        root = next(l for l in t.split('\n') if l.startswith('G * 0 '))
        f = root.split(' ')
        for field, tok in ((4, '0-30000000'), (6, '!+3-3000000000'),
                           (4, '0-18446744073709551615')):
            bad = _edit(t, root, ' '.join(f[:field] + [tok] + f[field + 1:]))
            tracemalloc.start()
            try:
                self.assert_bad(bad, 'beyond the 3 L records')
                peak = tracemalloc.get_traced_memory()[1]
            finally:
                tracemalloc.stop()
            self.assertLess(peak, 4 << 20, tok)
        for tok in ('V 0 0 AG 1,1 . . A:x:0-4000000000',):
            v = next(l for l in T.doc_text('linear').split('\n') if l.startswith('V '))
            self.assert_bad(_edit(T.doc_text('linear'), v, tok), 'beyond the')

    def test_names_keep_every_byte_but_lf(self):
        g = T.body('linear')
        odd = 'x\x0by\x0cz\x1c\x85   %0A\r\nend '
        g.labels[1].name = odd
        text = g.dump()
        self.assertEqual(text.count('\n'), int(text.split('\n')[-2][2:]))
        back = parse(text)
        self.assertEqual(odd, back.labels[1].name)
        self.assertEqual('acc3', back.labels[2].name)
        self.assertEqual(text, back.dump())

    def test_unknown_counter_round_trips(self):
        g = T.body('limits')
        self.assertEqual(7, g.arms['right'].counters['future_units'])
        self.assertEqual('future_units', list(g.arms['right'].counters)[-1])


if __name__ == '__main__':
    unittest.main()
