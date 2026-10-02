import json
import os
import random
import shlex
import socket
import subprocess
import time
import unittest

import requests

from base import PROTEIN_MODE, TestingBase, METAGRAPH

"""
End-to-end tests for the label-consistent traversal interface
(`metagraph traverse [--resolve]` and the server routes /resolve, /traverse,
/traverse/capabilities).

The fixture mirrors the format of the public `refseq33m` index: basic-mode
succinct graph at k=31, coordinate-aware annotation (row_diff_brwt_coord) with
one column per input file and a CoordToHeader (.seqs) sidecar, so labels are
per-sequence accession headers.

Sequence layout (random, fixed seed):
    acc1 = LEFT1 + ELEMENT + RIGHT1
    acc2 = LEFT1 + ELEMENT + RIGHT2      (same left flank, divergent right flank)
    acc3 = LEFT2 + ELEMENT + RIGHT1      (divergent left flank, same right flank)
all three in one file (one column), so the column label is the file and the
header labels are acc1/acc2/acc3.
"""

K = 31
BLOCK = 120


def _supports_traverse():
    res = subprocess.run([METAGRAPH, 'traverse'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = (res.stdout + res.stderr).decode()
    return 'Usage' in out and 'traverse' in out


@unittest.skipIf(PROTEIN_MODE, "traversal fixtures are DNA")
@unittest.skipUnless(_supports_traverse(), "`metagraph traverse` is not available in this build")
class TestTraverseBase(TestingBase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()

        rng = random.Random(20261002)

        def seq(n):
            return ''.join(rng.choice('ACGT') for _ in range(n))

        cls.left1, cls.left2 = seq(BLOCK), seq(BLOCK)
        cls.element = seq(BLOCK)
        cls.right1, cls.right2 = seq(BLOCK), seq(BLOCK)
        cls.records = {
            'acc1': cls.left1 + cls.element + cls.right1,
            'acc2': cls.left1 + cls.element + cls.right2,
            'acc3': cls.left2 + cls.element + cls.right1,
        }

        cls.fasta = cls.tempdir.name + '/ref.fa'
        with open(cls.fasta, 'w') as f:
            for name, s in cls.records.items():
                f.write(f'>{name}\n{s}\n')

        cls.graph = cls.tempdir.name + '/graph_k31.dbg'
        cls.anno_base = cls.tempdir.name + '/annotation'
        cls.anno = cls.anno_base + '.row_diff_brwt_coord.annodbg'

        cls._build_graph(cls.fasta, cls.graph, K, 'succinct', mode='basic')
        cls._annotate_graph(cls.fasta, cls.graph, cls.anno_base, 'row_diff_brwt_coord',
                            anno_type='filename')
        # Second pass: the CoordToHeader sidecar, so labels are the FASTA headers.
        # The loader strips the whole annotation extension and appends ".seqs", so the
        # sidecar must be written next to the annotation BASENAME, not next to the
        # fully-suffixed annotation file.
        cls._run_command(
            f'{METAGRAPH} annotate -i {cls.graph} --anno-filename --index-header-coords '
            f'-o {cls.anno_base} {cls.fasta}',
            'index-header-coords')
        assert os.path.exists(cls.anno_base + '.seqs'), 'CoordToHeader sidecar was not created'

    @classmethod
    def _traverse(cls, request, resolve=False):
        """Run the CLI on one request dict and return the parsed JSON."""
        path = cls.tempdir.name + '/request.json'
        with open(path, 'w') as f:
            json.dump(request, f)
        cmd = (f'{METAGRAPH} traverse {"--resolve" if resolve else ""} '
               f'-i {cls.graph} -a {cls.anno} {path}')
        res = subprocess.run(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out = res.stdout.decode()
        assert out.strip(), f'no output\nstderr: {res.stderr.decode()}'
        parsed = json.loads(out)
        if res.returncode and 'error' in parsed:
            print(f'  CLI error: {parsed["error"]}', flush=True)
        return parsed, res.returncode

    @staticmethod
    def _names(result):
        return {i: l['name'] for i, l in enumerate(result['label_dict'])}

    @staticmethod
    def _walks(arm, names):
        """leaf flank -> names of the labels at its end (a constrain-mode right arm).

        Only labels that travelled the spelled bases count (route_bp == 0); under
        `on_reconverge: keep` that is every label.
        """
        segments = {s['id']: s for s in arm['segments']}
        out = {}
        for path in arm['paths']:
            flank = ''.join(segments[sid]['sequence'] for sid in path['segments'])
            out[flank] = {names[e['label']] for e in path['end_labels'] if e['route_bp'] == 0}
        return out

    @staticmethod
    def _structural_walks(arm, names, permitted):
        """The maximal label-consistent walks of an ANNOTATE-mode right arm for the
        permitted set, read off the RECORDED label sets alone: the oracle side of the
        spec's §6.9 contract, sharing no label logic with the walker.
        """
        segments = {s['id']: s for s in arm['segments']}
        candidates = {}
        for path in arm['paths']:
            root = segments[path['segments'][0]]
            at = {0: {names[i] for i in root['labels']}}
            flank = ''
            for sid in path['segments']:
                seg = segments[sid]
                flank += seg['sequence']
                for run in seg['label_sets']:
                    assert not run['truncated'], 'a cut list makes the oracle unusable'
                    for d in range(run['from_bp'] + 1, run['to_bp'] + 1):
                        at[d] = {names[i] for i in run['labels']}
            assert len(at) == len(flank) + 1, 'the recorded runs do not cover the walk'
            for label in permitted:
                if label not in at[0]:
                    continue
                j = 0
                while j < len(flank) and label in at[j + 1]:
                    j += 1
                candidates.setdefault(label, set()).add(flank[:j])
        # maximal PER LABEL (a label whose walk ends inside another label's walk is a
        # claim of its own), then grouped by walk
        out = {}
        for label, walks in candidates.items():
            for w in walks:
                if not any(o != w and o.startswith(w) for o in walks):
                    out.setdefault(w, set()).add(label)
        return out


class TestTraverseCLI(TestTraverseBase):
    def test_resolve_header_labels(self):
        """resolve reports per-accession support and groups identical runs."""
        out, rc = self._traverse({
            'sequence': self.element,
            'discover': {'max_labels': 10, 'kind': 'header'},
        }, resolve=True)
        self.assertEqual(0, rc)
        self.assertEqual(K, out['k'])
        self.assertEqual('basic', out['regime'])
        self.assertEqual(BLOCK - K + 1, out['num_kmers'])
        names = sorted(l['label'] for l in out['labels'])
        self.assertEqual(['acc1', 'acc2', 'acc3'], names)
        for label in out['labels']:
            self.assertEqual('header', label['kind'])
            self.assertEqual(BLOCK - K + 1, label['kmers_supported'])
        # all three share the same support run -> one candidate with three labels
        self.assertEqual(1, len(out['candidates']))
        self.assertEqual(3, len(out['candidates'][0]['labels']))

    def test_resolve_select_freezes_seeds(self):
        out, rc = self._traverse({
            'sequence': self.left1 + self.element,
            'labels': ['acc1', 'acc2', 'acc3'],
            'select': {'policy': 'max_support', 'max_seeds': 2, 'min_block_bp': K},
        }, resolve=True)
        self.assertEqual(0, rc)
        seeds = out['selection']['seeds']
        self.assertTrue(seeds)
        # the whole query is carried by acc1 and acc2 (left1+element); acc3 only by element
        best = seeds[0]
        self.assertEqual(16, len(best['seed_id']))
        self.assertIn('acc1', best['labels'])
        self.assertIn('acc2', best['labels'])
        self.assertEqual(len(best['labels']), best['label_population']['included'])
        # the frozen seed is directly usable as a traverse seed
        self.assertTrue(best['sequence'])

    def test_traverse_splits_at_divergent_flank(self):
        """Seeding on the shared element: the right arm splits into RIGHT1/RIGHT2."""
        out, rc = self._traverse({
            'seeds': [{'sequence': self.element, 'labels': ['acc1', 'acc2', 'acc3']}],
            'strategy': {
                'direction': 'right',
                'bounds': {'max_extension_bp': BLOCK},
                'output': {'profile_bin_bp': 40},
            },
        })
        self.assertEqual(0, rc)
        result = out['results'][0]
        self.assertEqual(3, len(result['seed']['labels']))
        self.assertEqual([], result['seed']['dropped_labels'])
        right = result['arms']['right']
        self.assertEqual('complete', right['status'])
        self.assertNotIn('left', result['arms'])

        # two leaves: one carrying {acc1, acc3} (RIGHT1), one carrying {acc2} (RIGHT2)
        leaf_label_counts = sorted(p['n_labels'] for p in right['paths'])
        self.assertEqual([1, 2], leaf_label_counts)

        # reconstruct both flanks from the segments and compare to the sources
        flanks = set()
        segments = {s['id']: s for s in right['segments']}
        for path in right['paths']:
            flanks.add(''.join(segments[sid]['sequence'] for sid in path['segments']))
        self.assertEqual({self.right1[:BLOCK - 0][:len(self.right1)], self.right2},
                         {f[:len(self.right1)] for f in flanks})

        # a divergence is not an ambiguous branch: no label stopped with `branch`
        for seg in right['segments']:
            for ev in seg['events']:
                if ev['type'] == 'label_end':
                    self.assertNotEqual('branch', ev['reason'])

    def test_traverse_both_arms_and_label_summary(self):
        """acc3 diverges on the left, acc2 on the right; the summary shows both."""
        out, rc = self._traverse({
            'seeds': [{'sequence': self.element, 'labels': ['acc1', 'acc2', 'acc3']}],
            'strategy': {'bounds': {'max_extension_bp': BLOCK}},
        })
        self.assertEqual(0, rc)
        result = out['results'][0]
        self.assertIn('left', result['arms'])
        self.assertIn('right', result['arms'])
        names = [l['name'] for l in result['label_dict']]
        self.assertEqual(['acc1', 'acc2', 'acc3'], sorted(names))
        # every label reaches the full requested radius on at least one arm
        for entry in result['label_summary']:
            reach = max(entry['left']['reach_bp'], entry['right']['reach_bp'])
            self.assertGreater(reach, 0)

    def test_traverse_derives_labels_from_seed(self):
        """Omitting `labels` derives the permitted set from the seed (`permit: per_hit`)."""
        request = {
            'seeds': [{'sequence': self.element}],
            'strategy': {'direction': 'right', 'bounds': {'max_extension_bp': BLOCK},
                         'output': {'timing': False}},
        }
        out, rc = self._traverse(request)
        self.assertEqual(0, rc)
        # both new knobs are echoed in the normalized strategy
        labels = out['strategy']['labels']
        self.assertEqual(1000, labels['max_seed_labels'])
        self.assertEqual('auto', labels['seed_label_kind'])

        result = out['results'][0]
        seed = result['seed']
        self.assertTrue(seed['labels_from_seed'])
        self.assertEqual(3, seed['labels_supporting_total'])
        self.assertEqual(0, seed['labels_dropped'])
        self.assertEqual('', seed['labels_dropped_digest'])
        # the index carries a CoordToHeader, so the derived labels are the accessions
        self.assertEqual(['acc1', 'acc2', 'acc3'], sorted(seed['labels']))
        self.assertEqual(['header'] * 3, [l['kind'] for l in result['label_dict']])
        self.assertEqual([], seed['dropped_labels'])

        # The derived list is resubmittable VERBATIM, and naming it reproduces the whole
        # result object except `labels_from_seed` -- not just the arms. The derivation
        # reads WHOLE rows where an explicit list reads its own columns, so it maps more
        # coordinates; the number of rows it reads is the same one per seed k-mer.
        explicit = dict(request)
        explicit['seeds'] = [{'sequence': self.element, 'labels': seed['labels']}]
        named, rc = self._traverse(explicit)
        self.assertEqual(0, rc)
        self.assertFalse(named['results'][0]['seed']['labels_from_seed'])
        # an explicit list is its own support total, never 0
        self.assertEqual(3, named['results'][0]['seed']['labels_supporting_total'])

        def comparable(res):
            res = json.loads(json.dumps(res))
            res['seed'].pop('labels_from_seed')
            ann = res.pop('annotation')
            return res, ann

        mine, derived_ann = comparable(result)
        theirs, named_ann = comparable(named['results'][0])
        self.assertEqual(json.dumps(theirs, sort_keys=True), json.dumps(mine, sort_keys=True))
        self.assertEqual(named_ann['access_path'], derived_ann['access_path'])
        self.assertEqual(named_ann['rows_requested'], derived_ann['rows_requested'])
        self.assertEqual(named_ann['keys_mapped'], derived_ann['keys_mapped'])

        # the cap is reported, never silent
        request['strategy']['labels'] = {'max_seed_labels': 2}
        out, rc = self._traverse(request)
        self.assertEqual(0, rc)
        seed = out['results'][0]['seed']
        self.assertEqual(2, len(seed['labels']))
        self.assertEqual(3, seed['labels_supporting_total'])
        self.assertEqual(1, seed['labels_dropped'])
        self.assertEqual(16, len(seed['labels_dropped_digest']))

        # the echoed strategy must be resubmittable verbatim, "auto" included
        request['strategy']['labels'] = dict(labels)
        out, rc = self._traverse(request)
        self.assertEqual(0, rc, out.get('error'))
        self.assertEqual(labels, out['strategy']['labels'])

        # asking for the column kind derives the annotation column instead
        request['strategy']['labels'] = {'seed_label_kind': 'column'}
        out, rc = self._traverse(request)
        self.assertEqual(0, rc)
        self.assertEqual('column', out['strategy']['labels']['seed_label_kind'])
        result = out['results'][0]
        self.assertTrue(result['seed']['labels_from_seed'])
        self.assertEqual(1, len(result['seed']['labels']))
        self.assertEqual(['column'], [l['kind'] for l in result['label_dict']])

    def test_traverse_rejects_invalid_requests(self):
        for bad, expect in [
            ({'seeds': [{'sequence': self.element, 'labels': ['nope']}], 'strategy': {}},
             'nope'),
            # an omitted list derives; an explicitly empty one is an error naming it
            ({'seeds': [{'sequence': self.element, 'labels': []}], 'strategy': {}},
             'labels'),
            # a cost table is authored by name over the seed's own label list
            ({'seeds': [{'sequence': self.element}],
              'strategy': {'labels': {'change_cost': {'model': 'table', 'entries': []}}}},
             'change_cost'),
            ({'seeds': [{'sequence': 'ACGT', 'labels': ['acc1']}], 'strategy': {}},
             'k'),
            ({'seeds': [{'sequence': self.element, 'labels': ['acc1']}],
              'strategy': {'nonsense': 1}}, 'nonsense'),
            ({'seeds': [{'sequence': self.element, 'labels': ['acc1']}],
              'strategy': {'labels': {'extra': ['acc2']}}}, 'extra'),
            # the derived-set cap is bounded like every other knob
            ({'seeds': [{'sequence': self.element}],
              'strategy': {'labels': {'max_seed_labels': 100000000}}}, 'max_seed_labels'),
            ({'seeds': [{'sequence': self.element}],
              'strategy': {'labels': {'max_seed_labels': 0}}}, 'max_seed_labels'),
        ]:
            out, rc = self._traverse(bad)
            self.assertEqual(1, rc, f'expected failure for {bad}')
            self.assertIn('error', out)
            self.assertIn(expect, out['error'])

    def test_traverse_reports_caps(self):
        out, rc = self._traverse({
            'seeds': [{'sequence': self.element, 'labels': ['acc1']}],
            'strategy': {'direction': 'right', 'bounds': {'max_extension_bp': 10}},
        })
        self.assertEqual(0, rc)
        right = out['results'][0]['arms']['right']
        self.assertEqual('complete', right['status'])
        path = right['paths'][0]
        self.assertEqual(10, path['length_bp'])
        self.assertEqual('max_extension_bp', path['path_reason'])
        # a radius stop offers a continuation seed for iterative deepening
        self.assertIn('continuation', path)
        self.assertTrue(path['continuation']['sequence'])

    def test_traverse_derived_seed_failure_is_per_seed(self):
        """A seed no label carries in full fails ALONE, not the whole request.

        With explicit labels the caller picks the risk; with a derived set the failure is
        purely data-dependent and the caller has no lever, so one unlucky seed must not
        discard the other traversals of a batch.
        """
        # acc3 carries left2+element, acc2 carries element+right2, nobody carries both
        orphan = self.left2 + self.element + self.right2
        out, rc = self._traverse({
            'seeds': [{'seed_id': 'carried', 'sequence': self.element},
                      {'seed_id': 'orphan', 'sequence': orphan},
                      {'seed_id': 'carried2', 'sequence': self.left1 + self.element}],
            'strategy': {'direction': 'right', 'bounds': {'max_extension_bp': 10}},
        })
        self.assertEqual(0, rc, out.get('error'))
        self.assertEqual(3, len(out['results']))
        first, failed, third = out['results']
        for ok in (first, third):
            self.assertNotIn('error', ok)
            self.assertIn('right', ok['arms'])
        self.assertIn('error', failed)
        self.assertNotIn('arms', failed)
        self.assertEqual('orphan', failed['seed']['seed_id'])
        self.assertEqual(len(orphan), failed['seed']['length_bp'])
        self.assertTrue(failed['seed']['labels_from_seed'])
        self.assertIn('No label supports every k-mer', failed['error'])

        # an EXPLICIT label list is the caller's own choice, so it still fails the request
        out, rc = self._traverse({
            'seeds': [{'sequence': self.element, 'labels': ['acc1']},
                      {'sequence': orphan, 'labels': ['acc1', 'acc2', 'acc3']}],
            'strategy': {'direction': 'right', 'bounds': {'max_extension_bp': 10}},
        })
        self.assertEqual(1, rc)
        self.assertIn('error', out)

    def test_traverse_is_deterministic(self):
        request = {
            'seeds': [{'sequence': self.element, 'labels': ['acc1', 'acc2', 'acc3']}],
            'strategy': {'bounds': {'max_extension_bp': BLOCK},
                         'output': {'timing': False}},
        }
        first, _ = self._traverse(request)
        second, _ = self._traverse(request)
        self.assertEqual(json.dumps(first['results'], sort_keys=True),
                         json.dumps(second['results'], sort_keys=True))

    def test_traverse_annotate_mode_records_labels(self):
        """`labels.mode: annotate` follows every structural successor and records what is there.

        Seeding on the shared element, right arm: the structural trie forks into RIGHT1
        (carried by acc1 and acc3) and RIGHT2 (acc2). Nothing is filtered, the recorded
        sets name the carriers, and the constrained exhaustive run agrees leaf by leaf.
        """
        # a radius beyond the flanks, so that every walk ends at its real dead end
        radius = BLOCK + 10
        request = {
            'seeds': [{'sequence': self.element}],
            'strategy': {
                'exhaustive': True,
                'direction': 'right',
                'labels': {'mode': 'annotate', 'max_labels_per_node': 8},
                'bounds': {'max_extension_bp': radius},
                'output': {'timing': False},
            },
        }
        out, rc = self._traverse(request)
        self.assertEqual(0, rc, out.get('error'))
        self.assertIn('walk_rule', out)
        self.assertIn('edge twice', out['walk_rule'])
        self.assertIn('annotate', out['capabilities']['label_modes'])
        strategy = out['strategy']
        self.assertTrue(strategy['exhaustive'])
        self.assertEqual('annotate', strategy['labels']['mode'])
        self.assertEqual(8, strategy['labels']['max_labels_per_node'])
        self.assertEqual('unlimited', strategy['branching']['max_label_branches'])
        self.assertEqual('unlimited', strategy['branching']['max_splits_per_path'])
        self.assertEqual('keep', strategy['branching']['on_reconverge'])
        self.assertEqual('stop', strategy['frontier']['on_overflow'])

        result = out['results'][0]
        self.assertEqual('annotate', result['label_mode'])
        self.assertEqual([], result['seed']['labels'])
        names = {i: l['name'] for i, l in enumerate(result['label_dict'])}
        self.assertEqual({'acc1', 'acc2', 'acc3'}, set(names.values()))
        right = result['arms']['right']
        self.assertEqual('complete', right['status'])
        self.assertEqual(radius, right['complete_to_bp'])
        self.assertEqual({'cap': 8, 'max_seen': 3, 'nodes_truncated': 0}, right['labels_per_node'])
        # the trie view: one split right after the seed, one branch per flank, and the
        # branch counts (2 + 1) exceed nothing but need not sum to labels_before either
        self.assertEqual(1, len(right['splits']))
        split = right['splits'][0]
        self.assertEqual(0, split['prefix_bp'])
        self.assertEqual('divergence', split['kind'])
        self.assertEqual(3, split['labels_before'])
        self.assertEqual({1, 2}, {b['labels_distinct'] for b in split['branches']})
        for branch in split['branches']:
            self.assertFalse(branch['labels_truncated'])
            self.assertEqual(branch['labels_distinct'], len(branch['labels']))
        # every leaf ends structurally, carries no label ends, and its recorded sets are
        # exactly the carriers of that flank
        segments = {s['id']: s for s in right['segments']}
        recorded = {}
        for path in right['paths']:
            self.assertEqual('dead_end', path['path_reason'])
            self.assertEqual({}, path['end_reasons'])
            self.assertEqual([], path['end_labels'])
            flank = ''.join(segments[sid]['sequence'] for sid in path['segments'])
            alive = set(segments[path['segments'][0]]['labels'])
            for sid in path['segments']:
                for run in segments[sid]['label_sets']:
                    self.assertFalse(run['truncated'])
                    self.assertEqual(run['labels_total'], len(run['labels']))
                    alive &= set(run['labels'])
            recorded[flank] = {names[i] for i in alive}
        self.assertEqual({self.right1: {'acc1', 'acc3'}, self.right2: {'acc2'}}, recorded)

        # the echo is resubmittable verbatim and gives the same result (`output.timing`
        # is a request-level field the normalized strategy does not carry)
        request['strategy'] = dict(strategy)
        request['strategy'].pop('clamped')
        request['strategy']['output'] = dict(strategy['output'], timing=False)
        again, rc = self._traverse(request)
        self.assertEqual(0, rc, again.get('error'))
        self.assertEqual(json.dumps(out['results'], sort_keys=True),
                         json.dumps(again['results'], sort_keys=True))

        # the constrained exhaustive run over the carriers has the same leaves
        out, rc = self._traverse({
            'seeds': [{'sequence': self.element, 'labels': ['acc1', 'acc2', 'acc3']}],
            'strategy': {'exhaustive': True, 'direction': 'right',
                         'bounds': {'max_extension_bp': radius}, 'output': {'timing': False}},
        })
        self.assertEqual(0, rc, out.get('error'))
        result = out['results'][0]
        self.assertEqual('constrain', result['label_mode'])
        right = result['arms']['right']
        self.assertEqual('complete', right['status'])
        self.assertEqual(radius, right['complete_to_bp'])
        names = {i: l['name'] for i, l in enumerate(result['label_dict'])}
        segments = {s['id']: s for s in right['segments']}
        constrained = {}
        for path in right['paths']:
            flank = ''.join(segments[sid]['sequence'] for sid in path['segments'])
            constrained[flank] = {names[e['label']] for e in path['end_labels']}
        self.assertEqual(recorded, constrained)

    def test_traverse_exhaustive_rejects_conflicts_and_reports_the_complete_depth(self):
        """The preset refuses what would prune it; a tripped size cap reports complete_to_bp."""
        seed = {'sequence': self.element, 'labels': ['acc1', 'acc2', 'acc3']}
        for conflicting, field in [
                ({'branching': {'on_reconverge': 'merge'}}, 'on_reconverge'),
                ({'branching': {'max_label_branches': 2}}, 'max_label_branches'),
                ({'branching': {'max_splits_per_path': 3}}, 'max_splits_per_path'),
                ({'branching': {'min_successor_labels': 2}}, 'min_successor_labels'),
                ({'bounds': {'min_live_labels': 2}}, 'min_live_labels'),
                ({'frontier': {'on_overflow': 'beam'}}, 'on_overflow'),
                ({'labels': {'mode': 'annotate', 'extra': ['acc1']}}, 'labels.extra'),
                ({'labels': {'mode': 'annotate', 'loss_budget': 1}}, 'loss_budget'),
                ({'labels': {'mode': 'annotate'}, 'support': 'trace'}, 'support'),
                ({'labels': {'mode': 'annotate'}, 'branching': {'max_label_branches': 0}},
                 'max_label_branches'),
                ({'labels': {'mode': 'nonsense'}}, 'mode'),
                ({'branching': {'max_label_branches': 'lots'}}, 'max_label_branches'),
                # a pairwise cost keeps only the cheapest max_switch_sources sources,
                # which can leave a target label unentered: a number is refused
                ({'labels': {'change_cost': {'model': 'table', 'entries': [['acc1', 'acc2', 1]]},
                             'loss_budget': 1, 'max_switch_sources': 64}},
                 'max_switch_sources'),
        ]:
            strategy = {'exhaustive': True, 'direction': 'right'}
            strategy.update(conflicting)
            out, rc = self._traverse({'seeds': [seed], 'strategy': strategy})
            self.assertEqual(1, rc, f'accepted a conflict on {field}')
            self.assertIn(field, out['error'])
        # ... and "unlimited" is the preset's default for it, echoed as such
        out, rc = self._traverse({
            'seeds': [seed],
            'strategy': {'exhaustive': True, 'direction': 'right',
                         'labels': {'change_cost': {'model': 'table',
                                                    'entries': [['acc1', 'acc2', 1]]},
                                    'loss_budget': 1}},
        })
        self.assertEqual(0, rc, out.get('error'))
        self.assertEqual('unlimited', out['strategy']['labels']['max_switch_sources'])
        # a derived permitted set over max_seed_labels is refused per seed under the
        # preset (every walk of a dropped carrier would be missing), not cut
        out, rc = self._traverse({
            'seeds': [{'sequence': self.element}],
            'strategy': {'exhaustive': True, 'direction': 'right',
                         'labels': {'max_seed_labels': 2}},
        })
        self.assertEqual(0, rc)
        self.assertIn('max_seed_labels', out['results'][0]['error'])
        self.assertNotIn('arms', out['results'][0])
        # annotate mode has no permitted set: a seed label list is refused
        out, rc = self._traverse({'seeds': [seed],
                                  'strategy': {'labels': {'mode': 'annotate'}}})
        self.assertEqual(1, rc)
        self.assertIn('seeds[0].labels', out['error'])

        # a size cap trips mid-level: the partial level does not count, the status is
        # never "complete", and every walk up to complete_to_bp is present
        out, rc = self._traverse({
            'seeds': [{'sequence': self.element}],
            'strategy': {'exhaustive': True, 'direction': 'right',
                         'labels': {'mode': 'annotate'},
                         'bounds': {'max_extension_bp': BLOCK, 'max_steps': 15}},
        })
        self.assertEqual(0, rc, out.get('error'))
        # the annotate default of the per-node cap admits what the constrain side
        # admits as a derived set (max_seed_labels), so nothing is cut here ...
        self.assertEqual(1000, out['strategy']['labels']['max_labels_per_node'])
        right = out['results'][0]['arms']['right']
        self.assertEqual('truncated', right['status'])
        self.assertEqual('max_steps', right['cap_trigger']['reason'])
        # ... and the live-label counts at the trip are exact
        self.assertTrue(right['cap_trigger']['exact'])
        self.assertTrue(right['frontier_remaining']['exact'])
        self.assertTrue(all(g['exact'] for g in right['growth']))
        # the root segment carries the boundary node's true count (it is in no run)
        root = right['segments'][0]
        self.assertEqual(3, root['labels_total'])
        self.assertFalse(root['labels_truncated'])
        # two walks: 7 whole levels of two steps, the 8th level trips on its second head
        self.assertEqual(7, right['complete_to_bp'])
        segments = {s['id']: s for s in right['segments']}
        prefixes = set()
        for path in right['paths']:
            self.assertEqual('max_steps', path['path_reason'])
            flank = ''.join(segments[sid]['sequence'] for sid in path['segments'])
            self.assertGreaterEqual(len(flank), 7)
            prefixes.add(flank[:7])
        self.assertEqual({self.right1[:7], self.right2[:7]}, prefixes)

        # a cap that cuts: the counts over cut lists are lower bounds and say so
        out, rc = self._traverse({
            'seeds': [{'sequence': self.element}],
            'strategy': {'exhaustive': True, 'direction': 'right',
                         'labels': {'mode': 'annotate', 'max_labels_per_node': 1},
                         'bounds': {'max_extension_bp': BLOCK, 'max_steps': 15}},
        })
        self.assertEqual(0, rc, out.get('error'))
        right = out['results'][0]['arms']['right']
        self.assertGreater(right['labels_per_node']['nodes_truncated'], 0)
        self.assertFalse(right['cap_trigger']['exact'])
        self.assertFalse(right['frontier_remaining']['exact'])
        self.assertFalse(right['growth'][0]['exact'])
        root = right['segments'][0]
        self.assertEqual(3, root['labels_total'])
        self.assertEqual(1, len(root['labels']))
        self.assertTrue(root['labels_truncated'])

    def test_traverse_tuned_run_is_a_prefix_subset_of_the_exhaustive_trie(self):
        """A tuned run drops nothing silently (spec §6.9).

        Its leaves are prefixes of the exhaustive trie's under each of their labels, and
        the walk it omits carries a recorded reason: RIGHT2 is acc2 alone against acc1 and
        acc3 on RIGHT1, so a quorum of two prunes it with `minority` at the fork.
        """
        radius = BLOCK + 10
        seed = {'sequence': self.element, 'labels': ['acc1', 'acc2', 'acc3']}
        full, rc = self._traverse({
            'seeds': [seed],
            'strategy': {'exhaustive': True, 'direction': 'right',
                         'bounds': {'max_extension_bp': radius}},
        })
        self.assertEqual(0, rc, full.get('error'))
        names = self._names(full['results'][0])
        exhaustive = self._walks(full['results'][0]['arms']['right'], names)
        self.assertEqual({self.right1: {'acc1', 'acc3'}, self.right2: {'acc2'}}, exhaustive)

        tuned, rc = self._traverse({
            'seeds': [seed],
            'strategy': {'direction': 'right', 'bounds': {'max_extension_bp': radius},
                         'branching': {'min_successor_labels': 2, 'on_reconverge': 'keep'}},
        })
        self.assertEqual(0, rc, tuned.get('error'))
        right = tuned['results'][0]['arms']['right']
        self.assertEqual('complete', right['status'])
        names = self._names(tuned['results'][0])
        leaves = self._walks(right, names)
        # prefix-subset: every tuned leaf under each label extends to an exhaustive leaf
        for flank, labels in leaves.items():
            for label in labels:
                self.assertTrue(any(w.startswith(flank) and label in ls
                                    for w, ls in exhaustive.items()), (flank, label))
        # the omitted walk ...
        self.assertEqual({self.right1: {'acc1', 'acc3'}}, leaves)
        # ... has its reason on record at the fork: acc2 ended by the quorum, and the
        # branch event names both bases. Nothing split (one branch was followed), so
        # the root segment runs to the dead end and also holds acc1's and acc3's ends.
        root = [s for s in right['segments'] if not s['parents']][0]
        ends = [(ev['reason'], ev['at_bp'], names[ev['label']])
                for ev in root['events'] if ev['type'] == 'label_end']
        self.assertEqual([('minority', 0, 'acc2')], [e for e in ends if e[1] == 0])
        self.assertEqual({('dead_end', BLOCK, 'acc1'), ('dead_end', BLOCK, 'acc3')},
                         {e for e in ends if e[1] != 0})
        self.assertEqual(1, len(right['branch_events']))
        event = right['branch_events'][0]
        self.assertEqual(0, event['at_bp'])
        self.assertEqual({self.right1[0], self.right2[0]}, set(event['successors']))
        self.assertEqual([2, 1] if event['successors'][0] == self.right1[0] else [1, 2],
                         event['labels_per_successor'])
        self.assertEqual(['acc2'], [names[l] for l in event['dropped']])


class TestTraverseAPI(TestTraverseBase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.host = '127.0.0.1'
        os.environ['NO_PROXY'] = cls.host
        cls.port = 3601
        num_retries = 100
        while num_retries > 0:
            cls.server_process = cls._start_server(cls.graph, cls.anno, cls.host, cls.port)
            try:
                cls.server_process.wait(timeout=1)
            except subprocess.TimeoutExpired:
                break
            cls.port += 1
            num_retries -= 1
        if num_retries == 0:
            raise AssertionError("Couldn't start server")
        time.sleep(1)

    @classmethod
    def tearDownClass(cls):
        cls.server_process.kill()
        super().tearDownClass()

    @staticmethod
    def _start_server(graph, annotation, host, port):
        command = (f'{METAGRAPH} server_query -i {graph} -a {annotation} '
                   f'--port {port} --address {host} -p 2')
        return subprocess.Popen(shlex.split(command), shell=False,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    def _post(self, route, payload):
        return requests.post(url=f'http://{self.host}:{self.port}/{route}',
                             data=json.dumps(payload))

    def test_api_capabilities(self):
        ret = requests.get(url=f'http://{self.host}:{self.port}/traverse/capabilities')
        self.assertEqual(200, ret.status_code, ret.text)
        caps = ret.json()
        self.assertEqual(K, caps['k'])
        self.assertEqual('basic', caps['regime'])
        self.assertTrue(caps['has_coordinates'])
        self.assertTrue(caps['has_coord_to_header'])
        self.assertIn('forbid', caps['cost_models_available'])
        # A request may name no labels at all, which makes the SERVER choose the permitted
        # set: a deployment that configures nothing must not be the unlimited one.
        for cap in ('max_time_ms', 'max_seeds', 'max_seed_bp', 'max_seed_labels'):
            self.assertGreater(caps[cap], 0, cap)

    def test_api_enforces_server_caps(self):
        caps = requests.get(url=f'http://{self.host}:{self.port}/traverse/capabilities').json()
        strategy = {'direction': 'right', 'bounds': {'max_extension_bp': 10}}

        # a seed longer than the configured limit is refused before any graph work
        ret = self._post('traverse', {
            'seeds': [{'sequence': 'A' * (caps['max_seed_bp'] + 1)}], 'strategy': strategy})
        self.assertEqual(400, ret.status_code, ret.text)
        self.assertIn('longer than the server limit', ret.json()['error'])

        # ... and so are too many seeds in one request
        ret = self._post('traverse', {
            'seeds': [{'sequence': self.element}] * (caps['max_seeds'] + 1),
            'strategy': strategy})
        self.assertEqual(400, ret.status_code, ret.text)
        self.assertIn('seeds per request', ret.json()['error'])

        # the derived-set cap is clamped down, and the clamp is reported (never silent)
        over = caps['max_seed_labels'] + 1
        payload = dict(strategy)
        payload['labels'] = {'max_seed_labels': over}
        ret = self._post('traverse', {'seeds': [{'sequence': self.element}],
                                      'strategy': payload})
        self.assertEqual(200, ret.status_code, ret.text)
        out = ret.json()
        self.assertEqual(caps['max_seed_labels'], out['strategy']['labels']['max_seed_labels'])
        clamped = {c['field']: c for c in out['strategy']['clamped']}
        self.assertIn('labels.max_seed_labels', clamped)
        self.assertEqual(over, clamped['labels.max_seed_labels']['requested'])
        self.assertEqual(caps['max_seed_labels'],
                         clamped['labels.max_seed_labels']['effective'])

        # `time_budget_ms: 0` means "do not extend" for the WALK, but it is also the only
        # bound on deriving a permitted set the caller did not name, so a derived seed has
        # it clamped UP to the server cap rather than left unbounded.
        zero = {'direction': 'right',
                'bounds': {'max_extension_bp': 10, 'time_budget_ms': 0}}
        ret = self._post('traverse', {'seeds': [{'sequence': self.element}],
                                      'strategy': zero})
        self.assertEqual(200, ret.status_code, ret.text)
        out = ret.json()
        self.assertEqual(caps['max_time_ms'], out['strategy']['bounds']['time_budget_ms'])
        self.assertIn('bounds.time_budget_ms',
                      [c['field'] for c in out['strategy']['clamped']])
        # an explicit list keeps the zero: there it bounds nothing but the walk
        ret = self._post('traverse', {
            'seeds': [{'sequence': self.element, 'labels': ['acc1']}], 'strategy': zero})
        self.assertEqual(200, ret.status_code, ret.text)
        self.assertEqual(0, ret.json()['strategy']['bounds']['time_budget_ms'])
        self.assertEqual([], ret.json()['strategy']['clamped'])

    def test_api_resolve(self):
        ret = self._post('resolve', {
            'sequence': self.element,
            'labels': ['acc1', 'acc2', 'acc3'],
        })
        self.assertEqual(200, ret.status_code, ret.text)
        out = ret.json()
        self.assertEqual(3, len(out['labels']))
        self.assertEqual(1, len(out['candidates']))

    def test_api_traverse(self):
        ret = self._post('traverse', {
            'seeds': [{'sequence': self.element, 'labels': ['acc1', 'acc2', 'acc3']}],
            'strategy': {'direction': 'right', 'bounds': {'max_extension_bp': BLOCK}},
        })
        self.assertEqual(200, ret.status_code, ret.text)
        right = ret.json()['results'][0]['arms']['right']
        self.assertEqual('complete', right['status'])
        self.assertEqual(2, len(right['paths']))

    def test_api_traverse_without_labels(self):
        """A seed without `labels`: the permitted set comes from the seed itself."""
        ret = self._post('traverse', {
            'seeds': [{'sequence': self.element}],
            'strategy': {'direction': 'right', 'bounds': {'max_extension_bp': BLOCK}},
        })
        self.assertEqual(200, ret.status_code, ret.text)
        result = ret.json()['results'][0]
        self.assertTrue(result['seed']['labels_from_seed'])
        self.assertEqual(['acc1', 'acc2', 'acc3'], sorted(result['seed']['labels']))
        self.assertEqual(3, result['seed']['labels_supporting_total'])
        right = result['arms']['right']
        self.assertEqual('complete', right['status'])
        self.assertEqual(2, len(right['paths']))

    def test_api_invalid_request_is_400(self):
        ret = self._post('traverse', {'seeds': [], 'strategy': {}})
        self.assertEqual(400, ret.status_code)
        self.assertIn('seeds', ret.json()['error'])

    def test_api_routing_fields_are_accepted(self):
        """A 'graph' field must not trip the strict unknown-field check.

        The server consumes 'graph' / 'graph_path' to pick a shard before the request
        is parsed, so the parser has to declare them. It did not, which made every
        multi-graph request fail with "unknown field 'graph'" — invisible here until
        this test, because the rest of the suite runs a single-graph server.
        """
        for field in ('graph', 'graph_path'):
            ret = self._post('traverse', {
                field: 'anything',
                'seeds': [{'sequence': self.element, 'labels': ['acc1']}],
                'strategy': {'direction': 'right', 'bounds': {'max_extension_bp': 10}},
            })
            # single-graph server: naming a graph is a routing error, never a parse error
            self.assertEqual(400, ret.status_code, ret.text)
            self.assertNotIn('unknown field', ret.json()['error'],
                             f"'{field}' must be a known field")
            self.assertIn('single graph', ret.json()['error'])

    def test_api_rejects_access_override(self):
        """The access path follows the annotation; it is not client-selectable."""
        ret = self._post('traverse', {
            'seeds': [{'sequence': self.element, 'labels': ['acc1']}],
            'strategy': {'annotation': {'access': 'direct'}},
        })
        self.assertEqual(400, ret.status_code)
        self.assertIn('access', ret.json()['error'])

    def test_api_release_mismatch_is_400(self):
        ret = self._post('traverse', {
            'release': 'some-other-release',
            'seeds': [{'sequence': self.element, 'labels': ['acc1']}],
            'strategy': {},
        })
        # the server was started without --index-release, so a pinned request is
        # accepted only when the ids match; with no configured id it is ignored
        self.assertIn(ret.status_code, (200, 400))

    def _annotate(self, radius):
        ret = self._post('traverse', {
            'seeds': [{'sequence': self.element}],
            'strategy': {'exhaustive': True, 'direction': 'right',
                         'labels': {'mode': 'annotate'},
                         'bounds': {'max_extension_bp': radius}},
        })
        self.assertEqual(200, ret.status_code, ret.text)
        return ret.json()

    def test_api_traverse_annotate_mode(self):
        """HTTP, `labels.mode: annotate`: every structural successor is followed and the
        labels present are recorded; the response carries the completeness contract."""
        radius = BLOCK + 10
        out = self._annotate(radius)
        self.assertIn('edge twice', out['walk_rule'])
        self.assertEqual('traverse-0.2', out['algorithm_version'])
        self.assertTrue(out['strategy']['exhaustive'])
        self.assertEqual('annotate', out['strategy']['labels']['mode'])
        result = out['results'][0]
        self.assertEqual('annotate', result['label_mode'])
        self.assertEqual([], result['seed']['labels'])
        right = result['arms']['right']
        self.assertEqual('complete', right['status'])
        self.assertEqual(radius, right['complete_to_bp'])
        self.assertEqual(0, right['labels_per_node']['nodes_truncated'])
        self.assertEqual(3, right['labels_per_node']['max_seen'])
        for path in right['paths']:
            self.assertEqual('dead_end', path['path_reason'])
            self.assertEqual([], path['end_labels'])
            self.assertEqual({}, path['end_reasons'])
        names = self._names(result)
        structural = self._structural_walks(right, names, {'acc1', 'acc2', 'acc3'})
        self.assertEqual({self.right1: {'acc1', 'acc3'}, self.right2: {'acc2'}}, structural)
        # the trie view of the fork: one label per flank but acc1 and acc3 share RIGHT1
        split = right['splits'][0]
        self.assertEqual(0, split['prefix_bp'])
        self.assertEqual(3, split['labels_before'])
        self.assertEqual([1, 2], sorted(b['labels_distinct'] for b in split['branches']))
        # the preset refuses a knob that would prune it
        ret = self._post('traverse', {
            'seeds': [{'sequence': self.element}],
            'strategy': {'exhaustive': True, 'labels': {'mode': 'annotate'},
                         'frontier': {'on_overflow': 'beam'}},
        })
        self.assertEqual(400, ret.status_code, ret.text)
        self.assertIn('on_overflow', ret.json()['error'])

    def test_api_traverse_constrain_exhaustive_matches_the_structural_oracle(self):
        """HTTP, `labels.mode: constrain` with `exhaustive`: the label-constrained trie over a
        permitted set has exactly the structural oracle's leaves filtered by that set."""
        radius = BLOCK + 10
        structural = self._annotate(radius)['results'][0]
        oracle = self._structural_walks(structural['arms']['right'],
                                        self._names(structural), {'acc1', 'acc3'})
        # acc2's flank is not carried by the permitted set: it is not in E
        self.assertEqual({self.right1: {'acc1', 'acc3'}}, oracle)

        ret = self._post('traverse', {
            'seeds': [{'sequence': self.element, 'labels': ['acc1', 'acc3']}],
            'strategy': {'exhaustive': True, 'direction': 'right',
                         'bounds': {'max_extension_bp': radius}},
        })
        self.assertEqual(200, ret.status_code, ret.text)
        out = ret.json()
        self.assertEqual('constrain', out['strategy']['labels']['mode'])
        self.assertEqual('unlimited', out['strategy']['branching']['max_label_branches'])
        result = out['results'][0]
        self.assertEqual('constrain', result['label_mode'])
        right = result['arms']['right']
        self.assertEqual('complete', right['status'])
        self.assertEqual(radius, right['complete_to_bp'])
        self.assertEqual(oracle, self._walks(right, self._names(result)))
        # the path ends with the semantic reason the structural trie saw, a dead end,
        # reported per label in this mode (no path-level reason)
        self.assertEqual([{'dead_end': 2}], [p['end_reasons'] for p in right['paths']])
        self.assertNotIn('path_reason', right['paths'][0])


if __name__ == '__main__':
    unittest.main()
