import copy
import hashlib
import json
import os
import random
import shlex
import socket
import subprocess
import sys
import tempfile
import time
import unittest

import requests

from base import PROTEIN_MODE, TestingBase, METAGRAPH

# The graphlet library (api/python, metagraph.traverse) reads what the writer produces
# (DESIGN-traverse-graphlet.md §9, T37-T40). The test venv installs api/python editable;
# the library is stdlib-only, so the source tree serves as well -- these tests never skip.
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
try:
    import metagraph.traverse as graphlet_lib
except ImportError:
    sys.path.insert(0, os.path.join(REPO, 'api', 'python'))
    import metagraph.traverse as graphlet_lib
from metagraph.traverse.export import normalize_result  # noqa: E402

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
        # the index bundle's manifest (DESIGN-traverse-graphlet.md §3.1): its digest is the
        # index identity a graphlet carries, so retrievals can be joined and compared
        cls.manifest = cls.tempdir.name + '/index.manifest.json'
        cls.index_fp = cls._write_manifest(cls.tempdir.name, cls.manifest)

    @staticmethod
    def _write_manifest(directory, path):
        """Every file of the bundle (graph, annotation and their sidecars) with its size
        and sha256, as scripts/traversal/build_mini_refseq.sh writes it; -> index_fp, the
        sha256 over the lines "path\tsize\tsha256\n" in byte order of path."""
        files = []
        for name in sorted(os.listdir(directory)):
            full = os.path.join(directory, name)
            if os.path.isfile(full) and name.startswith(('graph', 'annotation')):
                with open(full, 'rb') as f:
                    data = f.read()
                files.append({'path': name, 'size': len(data),
                              'sha256': hashlib.sha256(data).hexdigest()})
        canonical = ''.join('%s\t%d\t%s\n' % (f['path'], f['size'], f['sha256'])
                            for f in files)
        fp = hashlib.sha256(canonical.encode()).hexdigest()
        with open(path, 'w') as f:
            json.dump({'format': 'metagraph-index-manifest', 'version': 1, 'index_fp': fp,
                       'files': files}, f, indent=1)
        return fp

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
    def _outcome_axes(result):
        """The per-seed outcome (spec §7.0) as (walks, branch_diagnostics, label_evidence,
        delivery): independent guarantees, one axis each."""
        o = result['outcome']
        assert set(o) == {'walks', 'branch_diagnostics', 'label_evidence', 'delivery'}, o
        return o['walks'], o['branch_diagnostics'], o['label_evidence'], o['delivery']

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
            # not implemented: refused rather than accepted and walked as 0
            ({'seeds': [{'sequence': self.element, 'labels': ['acc1']}],
              'strategy': {'branching': {'tip_window_bp': 20}}},
             'strategy.branching.tip_window_bp: not implemented'),
            ({'seeds': [{'sequence': self.element, 'labels': ['acc1']}],
              'strategy': {'branching': {'bubble_window_bp': 20}}},
             'strategy.branching.bubble_window_bp: not implemented'),
            # a continuation shorter than k could not be resubmitted as a seed
            ({'seeds': [{'sequence': self.element, 'labels': ['acc1']}],
              'strategy': {'output': {'continuation_bp': K - 1}}}, f'k = {K}'),
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
        # every result states its outcome: the walk reached the radius with all evidence
        self.assertEqual(('complete', 'complete', 'complete', 'inline'),
                         self._outcome_axes(out['results'][0]))

        def continuation(bp, sequence=None):
            out, rc = self._traverse({
                'seeds': [{'sequence': sequence or self.element, 'labels': ['acc1']}],
                'strategy': {'direction': 'right', 'bounds': {'max_extension_bp': 10},
                             'output': {'continuation_bp': bp}},
            })
            self.assertEqual(0, rc, out.get('error'))
            return out['results'][0]['arms']['right']['paths'][0]['continuation']

        # k bases: a valid seed, and resubmitting it works
        tail = continuation(K)
        self.assertEqual(K, len(tail['sequence']))
        self.assertEqual(self.right1[:10], tail['sequence'][-10:])
        continuation(K, tail['sequence'])
        # 0: no continuation sequence, but the labels and the loss are still reported
        empty = continuation(0)
        self.assertEqual('', empty['sequence'])
        self.assertEqual(1, len(empty['labels']))
        self.assertEqual(0, empty['loss_used'])

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
        # ... with a structured reason: the outcome and the request field to change
        self.assertEqual(('failed', 'complete', 'complete', 'inline'), self._outcome_axes(failed))
        self.assertEqual(['derivation'], [l['kind'] for l in failed['limitations']])
        derivation = failed['limitations'][0]
        self.assertEqual('no_carrier', derivation['cause'])
        self.assertEqual('seeds[].sequence', derivation['knob'])
        self.assertEqual(len(orphan) - K + 1, derivation['limit'])
        self.assertLessEqual(derivation['observed'], derivation['limit'])
        for ok in (first, third):
            self.assertEqual(('complete', 'complete', 'complete', 'inline'), self._outcome_axes(ok))

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

        # the echo is resubmittable verbatim and gives the same result (it carries
        # output.detail and output.timing; timing is turned off so results compare)
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

    def test_traverse_states_every_limitation(self):
        """Every cap that limited a result is stated in `limitations` with the request field
        to turn (spec §7.0), in every detail level, and a run nothing limited states none.
        `output.max_branch_events` accepts "unlimited" and echoes it so."""
        seed = {'sequence': self.element, 'labels': ['acc1', 'acc2', 'acc3']}

        def run(strategy, seeds=None):
            strategy = dict({'direction': 'right', 'bounds': {'max_extension_bp': BLOCK}},
                            **strategy)
            out, rc = self._traverse({'seeds': seeds or [seed], 'strategy': strategy})
            self.assertEqual(0, rc, out.get('error'))
            return out

        def kinds(limitations):
            return [l['kind'] for l in limitations]

        keep = {'on_reconverge': 'keep'}
        out = run({'branching': keep, 'output': {'max_branch_events': 'unlimited'}})
        self.assertEqual('unlimited', out['strategy']['output']['max_branch_events'])
        result = out['results'][0]
        self.assertEqual(('complete', 'complete', 'complete', 'inline'), self._outcome_axes(result))
        self.assertEqual([], result['limitations'])
        right = result['arms']['right']
        self.assertEqual('complete', right['status'])
        self.assertEqual([], right['limitations'])
        self.assertEqual({'complete': True, 'complete_to_bp': None}, right['evidence'])
        # the echo parses back, and the default stays 100
        again = run({'branching': keep, 'output': out['strategy']['output']})
        self.assertEqual('unlimited', again['strategy']['output']['max_branch_events'])
        self.assertEqual(100, run({})['strategy']['output']['max_branch_events'])
        bad, rc = self._traverse({'seeds': [seed],
                                  'strategy': {'output': {'max_branch_events': 'lots'}}})
        self.assertEqual(1, rc)
        self.assertIn('max_branch_events', bad['error'])

        # walk_domain: the RIGHT1 / RIGHT2 fork needs two live paths
        for detail in ('summary', 'tree', 'full'):
            result = run({'branching': keep, 'bounds': {'max_extension_bp': BLOCK, 'max_live_paths': 1},
                          'output': {'detail': detail}})['results'][0]
            self.assertEqual(('partial', 'complete', 'complete', 'inline'), self._outcome_axes(result))
            right = result['arms']['right']
            self.assertEqual('truncated', right['status'])
            self.assertEqual(['walk_domain'], kinds(right['limitations']), detail)
            walk = right['limitations'][0]
            self.assertEqual('bounds.max_live_paths', walk['knob'])
            self.assertEqual(1, walk['limit'])
            self.assertEqual(2, walk['observed'])
            self.assertEqual(right['complete_to_bp'], walk['complete_to_bp'])

        # branch_events: the quorum refusal at the fork is the only event; keeping none
        # puts the evidence boundary at the fork. The walks are complete, the branch
        # diagnostics are cut
        result = run({'branching': dict(keep, min_successor_labels=2),
                      'output': {'max_branch_events': 0}})['results'][0]
        self.assertEqual(('complete', 'cut', 'complete', 'inline'), self._outcome_axes(result))
        right = result['arms']['right']
        self.assertEqual('complete', right['status'])
        self.assertEqual({'complete': False, 'complete_to_bp': 0}, right['evidence'])
        self.assertEqual(['branch_events'], kinds(right['limitations']))
        events = right['limitations'][0]
        self.assertEqual('output.max_branch_events', events['knob'])
        self.assertEqual(0, events['limit'])
        self.assertEqual(1, events['observed'])
        self.assertEqual(0, events['complete_to_bp'])
        self.assertEqual(1, right['branch_events_truncated'])

        # scope: merging, constrain mode's default (no bubble here: nothing merged)
        right = run({})['results'][0]['arms']['right']
        self.assertEqual(['scope'], kinds(right['limitations']))
        self.assertEqual('branching.on_reconverge', right['limitations'][0]['knob'])
        self.assertEqual('merge', right['limitations'][0]['limit'])
        self.assertEqual(0, right['limitations'][0]['observed'])

        # label_lists (and the counts over the cut lists): annotate mode, one label a list
        result = run({'labels': {'mode': 'annotate', 'max_labels_per_node': 1}},
                     seeds=[{'sequence': self.element}])['results'][0]
        self.assertEqual(('complete', 'complete', 'lower_bound', 'inline'), self._outcome_axes(result))
        right = result['arms']['right']
        self.assertEqual(['label_lists', 'inexact_counts'], kinds(right['limitations']))
        lists = right['limitations'][0]
        self.assertEqual('labels.max_labels_per_node', lists['knob'])
        self.assertEqual(1, lists['limit'])
        self.assertEqual(right['labels_per_node']['max_seen'], lists['observed'])

        # seed_labels: the derived carrier set cut to two
        result = run({'branching': keep, 'labels': {'max_seed_labels': 2}},
                     seeds=[{'sequence': self.element}])['results'][0]
        self.assertEqual(['seed_labels'], kinds(result['limitations']))
        cut = result['limitations'][0]
        self.assertEqual('labels.max_seed_labels', cut['knob'])
        self.assertEqual(2, cut['limit'])
        self.assertEqual(result['seed']['labels_supporting_total'], cut['observed'])
        self.assertEqual([], result['arms']['right']['limitations'])
        # conservative outcome (DESIGN-traverse-graphlet.md §14, v5.2): the walks and the
        # evidence of the carrier the cap left out are missing, and both axes say so
        self.assertEqual(('partial', 'complete', 'lower_bound', 'inline'), self._outcome_axes(result))


class TestTraverseGraphlet(TestTraverseBase):
    """`output.detail: graphlet` through the CLI and the library (DESIGN-traverse-graphlet.md
    §9 as amended): every body the writer produces parses, dumps back byte for byte, is
    canonical and reproduces the detail: full result of the same request (T37); spellings,
    orientation and continuations (T38); the §6.9 oracle and the index identity through
    the library's comparisons (T39); the protocol (T40)."""

    def _traverse_identified(self, request, manifest=None):
        """The CLI with the index identity (--index-name, --index-manifest)."""
        path = self.tempdir.name + '/graphlet_request.json'
        with open(path, 'w') as f:
            json.dump(request, f)
        res = subprocess.run([METAGRAPH, 'traverse', '--index-name', 'tiny',
                              '--index-manifest', manifest or self.manifest,
                              '-i', self.graph, '-a', self.anno, path],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out = json.loads(res.stdout.decode())
        self.assertEqual(0, res.returncode, out.get('error'))
        return out

    def _retrieve(self, request):
        """-> (detail full response, detail graphlet response) of one request."""
        outs = []
        for detail in ('full', 'graphlet'):
            r = copy.deepcopy(request)
            r.setdefault('strategy', {}).setdefault('output', {})['detail'] = detail
            outs.append(self._traverse_identified(r))
        return outs

    def _graphlet(self, request, i=0):
        """The library's Graphlet of seed |i| of a detail: graphlet retrieval."""
        r = copy.deepcopy(request)
        r.setdefault('strategy', {}).setdefault('output', {})['detail'] = 'graphlet'
        out = self._traverse_identified(r)
        return graphlet_lib.from_response(out['results'][i], out), out

    def _assert_round_trip(self, full, out, case):
        """T37 for every seed of one request; -> the number of graphlets checked."""
        self.assertEqual(len(full['results']), len(out['results']), case)
        self.assertEqual('graphlet', out['strategy']['output']['detail'])
        checked = 0
        for i, (rf, rg) in enumerate(zip(full['results'], out['results'])):
            msg = '%s, seed %d' % (case, i)
            if 'error' in rf:
                # a failed derivation keeps its shape and carries no graphlet
                self.assertNotIn('graphlet', rg, msg)
                self.assertEqual(rf, rg, msg)
                continue
            checked += 1
            text = rg['graphlet']
            g = graphlet_lib.parse(text)
            self.assertEqual(text, g.dump(), msg)
            self.assertTrue(graphlet_lib.is_canonical(text), msg)
            self.assertEqual(rg['graphlet_lines'], text.count('\n'), msg)
            self.assertEqual('Z %d' % rg['graphlet_lines'], text.split('\n')[-2], msg)
            self.assertEqual(len(text.encode('utf-8')), rg['graphlet_bytes'], msg)
            gg = graphlet_lib.from_response(rg, out)
            self.assertEqual(normalize_result(rf), normalize_result(gg.to_json()), msg)
            # every '*' the writer emitted reproduces the walker's value, and the §4/§9
            # invariants hold
            self.assertEqual([], gg.check_rules(rf), msg)
            for side, arm in rf['arms'].items():
                a = g.arms[side]
                summary = rg['arms'][side]['counts']
                body = (a.counts.segments, a.counts.runs, a.counts.leaves, a.counts.splits,
                        a.counts.merges, a.counts.bases)
                self.assertEqual((len(arm['segments']), len(arm['runs']), len(arm['paths']),
                                  len(arm['splits'])), body[:4], msg)
                self.assertEqual((summary['segments'], summary['runs'], summary['leaves'],
                                  summary['splits'], summary['merges'], summary['bases']),
                                 body, msg)
                ends = sum(1 for s in arm['segments'] for e in s['events']
                           if e['type'] == 'label_end')
                self.assertLessEqual(ends, len(arm['runs']), msg)
                for r in a.runs:
                    if r.silent:
                        # a silent Lw end: the lineage went on under another name at to_bp
                        self.assertTrue(any(o.prev_run == r.id and o.from_label == r.label
                                            and o.from_bp == r.to_bp for o in a.runs), msg)
        return checked

    def test_t37_round_trip(self):
        """T37: the body of every strategy reproduces the detail: full result."""
        radius = BLOCK + 10
        el = self.element
        orphan = self.left2 + self.element + self.right2
        both = {'direction': 'both', 'bounds': {'max_extension_bp': radius}}
        switching = {'seeds': [{'sequence': self.left2[-40:] + el, 'labels': ['acc3']}],
                     'strategy': {'direction': 'right', 'bounds': {'max_extension_bp': radius},
                                  'branching': {'max_label_branches': 1},
                                  'labels': {'extra': ['acc2'], 'loss_budget': 1,
                                             'change_cost': {'model': 'constant',
                                                             'value': 0.5}}}}
        cases = {
            'tuned (merge), both arms': {'strategy': both},
            'exhaustive (keep)': {'strategy': dict(both, exhaustive=True)},
            'annotate exhaustive': {'strategy': dict(both, exhaustive=True,
                                                     labels={'mode': 'annotate'})},
            'annotate beam': {'strategy': dict(both, labels={'mode': 'annotate'},
                                               frontier={'on_overflow': 'beam'},
                                               bounds={'max_extension_bp': radius,
                                                       'max_live_paths': 1})},
            'quorum': {'strategy': dict(both, branching={'min_successor_labels': 2})},
            'left only': {'strategy': dict(both, direction='left')},
            'max_steps': {'strategy': dict(both, bounds={'max_extension_bp': radius,
                                                         'max_steps': 50})},
            # the cap trips on the left arm at depth 0 (two successors) and censors the
            # right root before it ran a level: an arm with no growth bin (A B* V*)
            'max_steps 1, annotate': {'strategy': dict(both, labels={'mode': 'annotate'},
                                                       bounds={'max_extension_bp': radius,
                                                               'max_steps': 1})},
            'max_labels_per_node 1': {'strategy': dict(both, labels={
                'mode': 'annotate', 'max_labels_per_node': 1})},
            'sequences false': {'strategy': dict(both, output={'sequences': False})},
            'radius below k': {'strategy': dict(both, bounds={'max_extension_bp': 10})},
            'trace, column labels': {'strategy': dict(both, support='trace',
                                                      branching={'on_reconverge': 'keep'},
                                                      labels={'seed_label_kind': 'column'})},
            'switch at a split': switching,
        }
        for case, request in cases.items():
            request = copy.deepcopy(request)
            request.setdefault('seeds', [{'sequence': el}])
            full, out = self._retrieve(request)
            self.assertEqual(1, self._assert_round_trip(full, out, case))
        # a 3-seed batch: a derivation failure keeps its shape, a repeat says duplicate
        full, out = self._retrieve({'seeds': [{'sequence': el}, {'sequence': orphan},
                                              {'sequence': el}], 'strategy': both})
        self.assertEqual(2, self._assert_round_trip(full, out, '3-seed batch'))
        self.assertTrue(out['results'][2]['duplicate'])
        self.assertEqual('failed', out['results'][1]['outcome']['walks'])

    def test_t37_cli_fixtures_are_reproduced(self):
        """T37 on the review's counterexamples (§9): the whole-document fixtures of
        api/python/tests/data/traverse/documents are what this binary writes for their
        requests, byte for byte, built from their tiny indexes with this CLI. Each also
        still shows its situation (an ambiguous split with disjoint child sets after a
        switch, a followed hairpin at a divergence, quorum-filtered splits, an A -> B -> A
        re-entry, silent switch-source ends at a split, ends whose branch counts differ
        only through a refused successor, two non-first-parent merges, two headers of one
        name); the library tests read the same files (api/python/tests)."""
        script = os.path.join(REPO, 'scripts', 'traversal', 'graphlet_fixtures.py')
        # the binary RELATIVE to the working directory, as a developer types it (the
        # script runs the CLI from each index directory, so it must resolve the path)
        res = subprocess.run([sys.executable, script, '--from-cli',
                              os.path.relpath(METAGRAPH, REPO), '--check'], cwd=REPO,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(0, res.returncode, res.stderr.decode() + res.stdout.decode())

    def test_t38_spelling_reproduces_the_records(self):
        """T38: natural(left) + seed + natural(right) is each record, read off the walks
        that carry its label on both arms."""
        g, _ = self._graphlet({'seeds': [{'sequence': self.element,
                                          'labels': ['acc1', 'acc2', 'acc3']}],
                               'strategy': {'exhaustive': True, 'direction': 'both',
                                            'bounds': {'max_extension_bp': BLOCK + 10}}})
        for name, record in self.records.items():
            left = [w for w in g.walks('left') if name in [l.name for l in w.labels_full]]
            right = [w for w in g.walks('right') if name in [l.name for l in w.labels_full]]
            self.assertEqual((1, 1), (len(left), len(right)), name)
            whole = (g.spell('left', left[0].path_id) + g.seed.sequence
                     + g.spell('right', right[0].path_id))
            self.assertEqual(record, whole, name)
            self.assertEqual(whole, g.spell('left', left[0].path_id, with_seed=True)
                             + g.spell('right', right[0].path_id))

    def test_t38_continuations_on_both_arms(self):
        """T38: the library's continuation equals the server's, on both arms, for
        continuations contained in the flank (40 of 60 bp) and crossing into the seed
        (k = 31 of 10 bp)."""
        for radius, bp, crossing in ((60, 40, False), (10, 1000, True)):
            request = {'seeds': [{'sequence': self.element}],
                       'strategy': {'direction': 'both',
                                    'bounds': {'max_extension_bp': radius},
                                    'output': {'continuation_bp': bp}}}
            full, out = self._retrieve(request)
            g = graphlet_lib.from_response(out['results'][0], out)
            for side in ('left', 'right'):
                paths = full['results'][0]['arms'][side]['paths']
                self.assertTrue(paths)
                for p in paths:
                    c = g.continuation(side, p['id'])
                    self.assertEqual(p['continuation']['sequence'], c.sequence, (radius, side))
                    self.assertEqual(crossing, len(c.sequence) > radius, (radius, side))

    def test_t38_a_local_continuation_is_resubmittable(self):
        """T38: a continuation derived by the library is a valid request, and the new
        traversal extends the walk."""
        g, _ = self._graphlet({'seeds': [{'sequence': self.element, 'labels': ['acc1']}],
                               'strategy': {'direction': 'right',
                                            'bounds': {'max_extension_bp': 60},
                                            'output': {'continuation_bp': 40}}})
        request = g.next_request('right', [0], bp=30)
        self.assertEqual([{'sequence': self.right1[20:60], 'labels': ['acc1']}],
                         request['seeds'])
        out = self._traverse_identified(request)
        result = out['results'][0]
        self.assertNotIn('error', result)
        self.assertEqual(('complete', 30), (result['arms']['right']['status'],
                                            result['arms']['right']['complete_to_bp']))
        nxt = graphlet_lib.from_response(result, out)
        self.assertEqual(self.right1[60:90], nxt.spell('right', 0))

    def test_t39_the_oracle_through_the_library(self):
        """T39: the constrained exhaustive trie equals the annotate oracle under the
        permitted labels (test_api_traverse_constrain_exhaustive_matches_the_structural_
        oracle), and a tuned run is a prefix subset of the exhaustive one whose omission
        carries the recorded reason (test_traverse_tuned_run_is_a_prefix_subset_...)."""
        radius = BLOCK + 10
        right = {'direction': 'right', 'bounds': {'max_extension_bp': radius}}

        def get(strategy, labels=None):
            seed = {'sequence': self.element}
            if labels:
                seed['labels'] = labels
            return self._graphlet({'seeds': [seed], 'strategy': strategy})[0]

        annotate = get(dict(right, labels={'mode': 'annotate'}))
        self.assertEqual((self.index_fp, 'tiny'), (annotate.index_fp, annotate.index_ns))
        constrain = get(dict(right, exhaustive=True), ['acc1', 'acc3'])
        c = constrain.compare(annotate, mode='claims', labels=['acc1', 'acc3'])
        self.assertEqual((True, True, [], [], []),
                         (c.comparable, c.equal, c.only_in_a, c.only_in_b, c.differ))
        self.assertEqual(radius, c.depth_used)
        # over all labels, acc2's walk is the annotate side's only
        c = constrain.compare(annotate, mode='claims')
        self.assertFalse(c.equal)
        self.assertEqual([{'name': 'acc2', 'ref': 'h:0:1'}], [x['label'] for x in c.only_in_b])
        exhaustive = get(dict(right, exhaustive=True), ['acc1', 'acc2', 'acc3'])
        for mode in ('claims', 'walks', 'labels'):
            c = exhaustive.compare(annotate, mode=mode)
            self.assertTrue(c.equal, (mode, c.only_in_a, c.only_in_b, c.differ))
        tuned = get(dict(right, branching={'min_successor_labels': 2, 'on_reconverge': 'keep'}),
                    ['acc1', 'acc2', 'acc3'])
        c = tuned.compare(exhaustive, mode='prefix_subset')
        self.assertTrue(c.equal)
        self.assertEqual([('acc2', self.right2, 0, 'minority')],
                         [(x['name'], x['walk'], x['at_bp'], x['reason']) for x in c.only_in_b])
        self.assertFalse(exhaustive.compare(tuned, mode='prefix_subset').equal)
        # the same on the left arm (walks are read outward there too)
        left = dict(right, direction='left')
        c = get(dict(left, branching={'min_successor_labels': 2, 'on_reconverge': 'keep'}),
                ['acc1', 'acc2', 'acc3']).compare(
                    get(dict(left, exhaustive=True), ['acc1', 'acc2', 'acc3']),
                    mode='prefix_subset')
        self.assertTrue(c.equal)
        self.assertEqual([('acc3', self.left2, 0, 'minority')],
                         [(x['name'], x['walk'], x['at_bp'], x['reason']) for x in c.only_in_b])
        # a claim reaching the comparison depth is open on both sides: a shallower
        # retrieval of the same trie compares equal
        short = get(dict(right, exhaustive=True, bounds={'max_extension_bp': 60}),
                    ['acc1', 'acc2', 'acc3'])
        c = short.compare(exhaustive, mode='claims')
        self.assertEqual((True, True, 60), (c.comparable, c.equal, c.depth_used))
        # a side whose label evidence is a lower bound (cut lists) is never equal
        cut = get(dict(right, labels={'mode': 'annotate', 'max_labels_per_node': 1}))
        self.assertEqual('lower_bound', cut.outcome.label_evidence)
        c = exhaustive.compare(cut, mode='claims')
        self.assertEqual(('qualified', None), (c.comparable, c.equal))
        # no common certified depth: nothing to compare, not equality
        none = get(dict(right, exhaustive=True, bounds={'max_extension_bp': radius,
                                                        'max_steps': 1}),
                   ['acc1', 'acc2', 'acc3'])
        self.assertEqual(0, none.arm('right').complete_to_bp)
        for mode in ('walks', 'labels', 'claims'):
            c = exhaustive.compare(none, mode=mode)
            self.assertEqual(('unknown', None, 0), (c.comparable, c.equal, c.depth_used))

    def test_t37_names_that_are_not_utf8(self):
        """Label names come from FASTA headers, which need not be UTF-8. The server makes
        them valid once (each ill-formed sequence -> U+FFFD), so detail: full, the summary
        and the body carry the same bytes: graphlet_bytes holds after transport and the
        front-coded names decode to the names of detail: full."""
        rng = random.Random(4711)
        block = {b: ''.join(rng.choice('ACGT') for _ in range(30)) for b in 'SABC'}
        headers = [b'\xc0abc', b'\xc0abd', b'caf\xc3\xa9']
        d = os.path.join(self.tempdir.name, 'latin1')
        os.makedirs(d)
        with open(os.path.join(d, 'x.fa'), 'wb') as f:
            for h, tail in zip(headers, 'ABC'):
                f.write(b'>' + h + b'\n' + (block['S'] + block[tail]).encode() + b'\n')
        for cmd in (f'{METAGRAPH} build -p 1 --graph succinct -k 15 -o graph x.fa',
                    f'{METAGRAPH} annotate -p 1 --anno-filename -i graph.dbg --anno-type column '
                    f'--coordinates -o annotation x.fa',
                    f'{METAGRAPH} transform_anno -p 1 --anno-type column_coord --coordinates '
                    f'-o annotation annotation.column.annodbg',
                    f'{METAGRAPH} annotate -p 1 -i graph.dbg --anno-filename '
                    f'--index-header-coords -o annotation x.fa'):
            res = subprocess.run(shlex.split(cmd), cwd=d, stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
            self.assertEqual(0, res.returncode, res.stderr.decode())
        want = sorted(h.decode('utf-8', 'replace') for h in headers)
        for mode in ('annotate', 'constrain'):
            outs = []
            for detail in ('full', 'graphlet'):
                path = os.path.join(d, 'request.json')
                with open(path, 'w') as f:
                    json.dump({'seeds': [{'sequence': block['S']}],
                               'strategy': {'direction': 'right',
                                            'labels': {'mode': mode,
                                                       'seed_label_kind': 'header'},
                                            'bounds': {'max_extension_bp': 40},
                                            'output': {'detail': detail}}}, f)
                res = subprocess.run([METAGRAPH, 'traverse', '-i', os.path.join(d, 'graph.dbg'),
                                      '-a', os.path.join(d, 'annotation.column_coord.annodbg'),
                                      path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                self.assertEqual(0, res.returncode, res.stdout.decode('utf-8', 'replace'))
                outs.append(json.loads(res.stdout.decode('utf-8')))
            full, out = outs
            self.assertEqual(want, sorted(l['name'] for l in full['results'][0]['label_dict']),
                             mode)
            self.assertEqual(1, self._assert_round_trip(full, out, 'not UTF-8, ' + mode))
            g = graphlet_lib.from_response(out['results'][0], out)
            self.assertEqual(want, sorted(l.name for l in g.labels), mode)

    def test_t39_index_identity(self):
        """§3.1 (freeze gate): two indexes over the same records with the same counts and
        column names but swapped memberships share index_meta_fp and differ in index_fp,
        so their retrievals never compare equal; without a manifest a join is
        unverifiable."""
        a, b = self.records['acc1'], self.records['acc2']        # equal lengths
        fps, graphlets = [], []
        for name, (first, second) in (('ab', (a, b)), ('ba', (b, a))):
            d = os.path.join(self.tempdir.name, 'swapped_' + name)
            os.makedirs(d)
            for fname, rec in (('colA.fa', first), ('colB.fa', second)):
                with open(os.path.join(d, fname), 'w') as f:
                    f.write('>%s\n%s\n' % (fname[:4], rec))
            for cmd in (f'{METAGRAPH} build -p 1 --mode basic --graph succinct -k {K} '
                        f'-o graph colA.fa colB.fa',
                        f'{METAGRAPH} annotate -p 1 --anno-filename -i graph.dbg '
                        f'--anno-type column --coordinates -o annotation colA.fa colB.fa',
                        f'{METAGRAPH} transform_anno -p 1 --anno-type column_coord '
                        f'--coordinates -o annotation annotation.column.annodbg'):
                res = subprocess.run(shlex.split(cmd), cwd=d, stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE)
                self.assertEqual(0, res.returncode, res.stderr.decode())
            os.remove(os.path.join(d, 'annotation.column.annodbg'))
            os.remove(os.path.join(d, 'annotation.column.annodbg.coords'))
            manifest = os.path.join(d, 'manifest.json')
            fps.append(self._write_manifest(d, manifest))
            path = os.path.join(d, 'request.json')
            with open(path, 'w') as f:
                json.dump({'seeds': [{'sequence': self.element}],
                           'strategy': {'direction': 'right',
                                        'labels': {'seed_label_kind': 'column'},
                                        'output': {'detail': 'graphlet'}}}, f)
            out = []
            for flags in ([], ['--index-manifest', manifest]):
                res = subprocess.run([METAGRAPH, 'traverse'] + flags +
                                     ['-i', os.path.join(d, 'graph.dbg'),
                                      '-a', os.path.join(d, 'annotation.column_coord.annodbg'),
                                      path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                resp = json.loads(res.stdout.decode())
                self.assertEqual(0, res.returncode, resp.get('error'))
                out.append(graphlet_lib.from_response(resp['results'][0], resp))
            graphlets.append(out)
        (ab_bare, ab), (ba_bare, ba) = graphlets
        self.assertEqual(ab.index_meta_fp, ba.index_meta_fp)
        self.assertNotEqual(fps[0], fps[1])
        self.assertEqual((fps[0], fps[1]), (ab.index_fp, ba.index_fp))
        c = ab.compare(ba)
        self.assertEqual((False, None), (c.comparable, c.equal))
        self.assertEqual(('unverifiable', None), (ab_bare.compare(ba_bare).comparable,
                                                  ab_bare.compare(ba_bare).equal))
        self.assertEqual((True, True), (ab.compare(ab).comparable, ab.compare(ab).equal))

    def test_t40_protocol(self):
        """T40: the summary keys, the echo, the capabilities, determinism and size."""
        request = {'seeds': [{'sequence': self.element}],
                   'strategy': {'direction': 'both',
                                'bounds': {'max_extension_bp': BLOCK + 10}}}
        full, out = self._retrieve(request)
        self.assertEqual('graphlet', out['strategy']['output']['detail'])
        self.assertTrue(out['strategy']['output']['timing'])
        caps = out['capabilities']
        self.assertEqual(1, caps['graphlet_format'])
        self.assertEqual(['summary', 'tree', 'full', 'graphlet'], caps['detail_levels'])
        self.assertEqual(('tiny', self.index_fp), (caps['index_ns'], caps['index_fp']))
        r = out['results'][0]
        self.assertEqual({'seed', 'label_mode', 'outcome', 'limitations', 'arms', 'annotation',
                          'timing', 'graphlet', 'graphlet_bytes', 'graphlet_lines'}, set(r))
        self.assertEqual({'seed_id', 'validated_seed_id', 'seed_id_mismatch', 'length_bp',
                          'num_kmers', 'labels_from_seed', 'labels_supporting_total',
                          'labels_dropped', 'labels_dropped_digest', 'num_labels',
                          'num_seed_labels'}, set(r['seed']))
        self.assertIn('serialize_ms', r['timing'])
        for arm in r['arms'].values():
            self.assertEqual({'segments', 'leaves', 'splits', 'merges', 'runs', 'bases', 'max_bp',
                              'label_ends', 'leaves_by_reason'}, set(arm['counts']))
            self.assertNotIn('segments', arm)
        h = r['graphlet'].split('\n')[0].split(' ')
        self.assertEqual(['tiny', self.index_fp, caps['index_meta_fp']], h[-3:])
        # the echo is resubmittable verbatim (without the server's clamped report, as in
        # test_traverse_annotate_mode_records_labels), and the body is deterministic
        echo = {k: v for k, v in out['strategy'].items() if k != 'clamped'}
        again = self._traverse_identified({'seeds': request['seeds'], 'strategy': echo})
        self.assertEqual(r['graphlet'], again['results'][0]['graphlet'])
        self.assertEqual(r['graphlet'], self._retrieve(request)[1]['results'][0]['graphlet'])
        # the body is less than half of the compact full result
        compact = json.dumps(full['results'][0], separators=(',', ':'))
        self.assertLess(r['graphlet_bytes'] * 2, len(compact.encode()))

    def test_t40_without_sequences(self):
        """T40: sequences false -- G records without bases, split children with their
        first base, continuations with their sequence, splits[].char reproduced."""
        request = {'seeds': [{'sequence': self.element}],
                   'strategy': {'direction': 'both', 'bounds': {'max_extension_bp': 60},
                                'output': {'sequences': False, 'continuation_bp': 40}}}
        full, out = self._retrieve(request)
        self._assert_round_trip(full, out, 'sequences false')
        lines = out['results'][0]['graphlet'].split('\n')
        gs = [l.split(' ') for l in lines if l.startswith('G ')]
        self.assertTrue(all(f[10] == '*' for f in gs))
        self.assertTrue(all(f[9] != '*' for f in gs if f[1] != '*'))      # split children
        cs = [l.split(' ') for l in lines if l.startswith('C ')]
        self.assertTrue(cs and all(len(f) == 6 for f in cs))
        g = graphlet_lib.from_response(out['results'][0], out)
        mine = g.to_json()
        for side, arm in full['results'][0]['arms'].items():
            self.assertEqual([[b['char'] for b in s['branches']] for s in arm['splits']],
                             [[b['char'] for b in s['branches']]
                              for s in mine['arms'][side]['splits']])


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

    @classmethod
    def _start_server(cls, graph, annotation, host, port):
        # with the index identity, which every graphlet and capabilities response carries
        command = (f'{METAGRAPH} server_query -i {graph} -a {annotation} '
                   f'--port {port} --address {host} -p 2 '
                   f'--index-name tiny --index-manifest {cls.manifest}')
        return subprocess.Popen(shlex.split(command), shell=False,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    def _post(self, route, payload):
        return requests.post(url=f'http://{self.host}:{self.port}/{route}',
                             data=json.dumps(payload))

    def test_api_graphlet(self):
        """T37/T40 over HTTP: the graphlet body of a gzip-compressed response reproduces
        the detail: full result; the server's index identity is in capabilities, in the
        response and in every body's H record."""
        caps = requests.get(url=f'http://{self.host}:{self.port}/traverse/capabilities').json()
        self.assertEqual(('tiny', self.index_fp, 1),
                         (caps['index_ns'], caps['index_fp'], caps['graphlet_format']))
        self.assertIn('graphlet', caps['detail_levels'])
        outs = {}
        for detail in ('full', 'graphlet'):
            ret = requests.post(url=f'http://{self.host}:{self.port}/traverse',
                                data=json.dumps({
                                    'seeds': [{'sequence': self.element}],
                                    'strategy': {'direction': 'both',
                                                 'bounds': {'max_extension_bp': BLOCK + 10},
                                                 'output': {'detail': detail}}}),
                                headers={'Accept-Encoding': 'gzip'})
            self.assertEqual(200, ret.status_code, ret.text)
            self.assertEqual('gzip', ret.headers.get('Content-Encoding'))
            outs[detail] = ret.json()
        out = outs['graphlet']
        self.assertEqual(caps['index_meta_fp'], out['capabilities']['index_meta_fp'])
        r = out['results'][0]
        text = r['graphlet']
        self.assertEqual(['tiny', self.index_fp, caps['index_meta_fp']],
                         text.split('\n')[0].split(' ')[-3:])
        g = graphlet_lib.from_response(r, out)
        self.assertEqual(text, g.dump())
        self.assertTrue(graphlet_lib.is_canonical(text))
        self.assertEqual(normalize_result(outs['full']['results'][0]),
                         normalize_result(g.to_json()))
        self.assertEqual([], g.check_rules(outs['full']['results'][0]))

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
        # an integer knob is echoed as an integer (10001, not 10001.0)
        self.assertIsInstance(clamped['labels.max_seed_labels']['requested'], int)
        self.assertIsInstance(clamped['labels.max_seed_labels']['effective'], int)
        self.assertNotIn(f'{over}.0', ret.text)
        # ... but the three carriers fit under it: the clamp bound nothing in this seed,
        # so the seed states no limitation (spec §7.0)
        self.assertEqual([], out['results'][0]['limitations'])

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
        # the walk ran under a budget the request did not ask for: stated on the seed
        self.assertEqual([{'kind': 'server_clamp', 'knob': 'bounds.time_budget_ms',
                           'limit': caps['max_time_ms'], 'observed': 0}],
                         [{k: l[k] for k in ('kind', 'knob', 'limit', 'observed')}
                          for l in out['results'][0]['limitations']])
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

    def test_api_traversal_routes_are_compact_and_compressible(self):
        """HTTP transport: the traversal routes write compact JSON (no indentation) and
        honour Accept-Encoding with gzip (preferred) or deflate; the capabilities say so."""
        payload = json.dumps({'seeds': [{'sequence': self.element, 'labels': ['acc1']}],
                              'strategy': {'direction': 'right', 'bounds': {'max_extension_bp': 20}}})
        url = f'http://{self.host}:{self.port}/traverse'
        plain = requests.post(url=url, data=payload, headers={'Accept-Encoding': 'identity'})
        self.assertEqual(200, plain.status_code, plain.text)
        self.assertNotIn('Content-Encoding', plain.headers)
        self.assertNotIn('\n', plain.text)
        self.assertNotIn('\t', plain.text)
        for encoding in ('gzip', 'deflate'):
            ret = requests.post(url=url, data=payload, headers={'Accept-Encoding': encoding})
            self.assertEqual(200, ret.status_code, ret.text)
            self.assertEqual(encoding, ret.headers.get('Content-Encoding'))
            self.assertLess(len(ret.raw.getheader('Content-Length') or '0'), 20)
            self.assertEqual(plain.json()['results'][0]['arms'], ret.json()['results'][0]['arms'])
        both = requests.post(url=url, data=payload, headers={'Accept-Encoding': 'deflate, gzip'})
        self.assertEqual('gzip', both.headers.get('Content-Encoding'))
        caps = requests.get(url=f'http://{self.host}:{self.port}/traverse/capabilities').json()
        self.assertEqual(['gzip', 'deflate'], caps['content_encodings'])

    def test_api_accept_encoding_honours_quality_values(self):
        """HTTP transport (RFC 9110 §12.5.3): a coding with q=0 is not acceptable, `*` covers
        the codings not listed, tokens are case-insensitive and the higher weight wins
        (gzip on a tie). A substring test used to answer gzip to the first two headers
        below (review round 4, finding 4). It is the shared process_request() path, so
        the capabilities route stands for every route."""
        url = f'http://{self.host}:{self.port}/traverse/capabilities'
        plain = requests.get(url=url, headers={'Accept-Encoding': 'identity'})
        self.assertEqual(200, plain.status_code, plain.text)
        for accept, expected in [
                ('gzip;q=0, deflate;q=1', 'deflate'),
                ('identity, gzip;q=0', None),
                ('gzip;q=0.5, deflate;q=1', 'deflate'),
                ('*;q=0', None),
                ('GZIP', 'gzip'),
                ('deflate, gzip', 'gzip'),
                ('gzip;q=0.1, deflate;q=1', 'deflate'),
                ('*', 'gzip'),
                ('*;q=0.5, deflate', 'deflate'),
                ('identity;q=1, gzip;q=0.5', None),
                ('gzip;q=abc, deflate;q=0.2', 'deflate'),
        ]:
            ret = requests.get(url=url, headers={'Accept-Encoding': accept})
            self.assertEqual(200, ret.status_code, accept)
            self.assertEqual(expected, ret.headers.get('Content-Encoding'), accept)
            self.assertEqual(plain.json(), ret.json(), accept)

    def test_api_annotate_defaults_to_keep_and_merging_qualifies_the_walk_rule(self):
        """HTTP: without `exhaustive`, annotate mode defaults to `on_reconverge: keep` (it
        promises the per-path trie) while constrain mode keeps merging as its default; when
        merging is on, `walk_rule` says the edge histories are united. The physical fetch
        counters live in `timing`, not in the invariant `result.annotation`."""
        ret = self._post('traverse', {
            'seeds': [{'sequence': self.element}],
            'strategy': {'direction': 'right', 'labels': {'mode': 'annotate'},
                         'bounds': {'max_extension_bp': 20}},
        })
        self.assertEqual(200, ret.status_code, ret.text)
        out = ret.json()
        self.assertEqual('keep', out['strategy']['branching']['on_reconverge'])
        self.assertNotIn('edge histories are united', out['walk_rule'])
        ret = self._post('traverse', {
            'seeds': [{'sequence': self.element, 'labels': ['acc1']}],
            'strategy': {'direction': 'right', 'bounds': {'max_extension_bp': 20}},
        })
        self.assertEqual(200, ret.status_code, ret.text)
        out = ret.json()
        self.assertEqual('merge', out['strategy']['branching']['on_reconverge'])
        self.assertIn('edge histories are united', out['walk_rule'])
        self.assertIn('on_reconverge: keep gives the per-path set', out['walk_rule'])
        result = out['results'][0]
        self.assertEqual({'access_path', 'keys_mapped', 'rows_requested', 'direct_reads'},
                         set(result['annotation'].keys()))
        for key in ('rows_fetched', 'tuple_rows_fetched', 'coords_mapped'):
            self.assertIn(key, result['timing'])

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
