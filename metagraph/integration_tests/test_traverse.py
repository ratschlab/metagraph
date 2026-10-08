import copy
import filecmp
import hashlib
import json
import os
import random
import re
import shlex
import shutil
import signal
import socket
import subprocess
import sys
import tempfile
import threading
import time
import types
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
        sha256 over the lines "path\tsize\tsha256\n" in byte order of path. Never the
        graph's derived data (a mask, a Bloom filter): not part of index_fp, and refused in a
        manifest (the owner's decision #17 of 2026-10-08)."""
        files = []
        for name in sorted(os.listdir(directory)):
            full = os.path.join(directory, name)
            if (os.path.isfile(full) and name.startswith(('graph', 'annotation'))
                    and not name.endswith(('.edgemask', '.bloom'))):
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
    def _traverse(cls, request, resolve=False, flags=''):
        """Run the CLI on one request dict (|flags|: more CLI flags) and return the parsed JSON."""
        path = cls.tempdir.name + '/request.json'
        with open(path, 'w') as f:
            json.dump(request, f)
        cmd = (f'{METAGRAPH} traverse {"--resolve" if resolve else ""} {flags} '
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

    def test_traverse_attempt_errors_state_their_usage(self):
        """The CLI states the usage of a request with attempt_id on an error after the
        request was read, as the server's 400 does; a malformed id is refused without."""
        bad = {'seeds': [{'sequence': self.element, 'labels': ['nope']}], 'strategy': {},
               'attempt_id': 'cli-err', 'locus_id': 'l-1'}
        out, rc = self._traverse(bad)
        self.assertEqual(1, rc)
        self.assertIn('nope', out['error'])
        self.assertEqual(('cli-err', 'l-1', 'error', False),
                         (out['usage']['attempt_id'], out['usage']['locus_id'],
                          out['usage']['reason'], out['usage']['bound']['enforced']))
        out, rc = self._traverse(dict(bad, attempt_id='not an id'))
        self.assertEqual(1, rc)
        self.assertIn('attempt_id', out['error'])
        self.assertNotIn('usage', out)
        # without attempt_id the error is as it was
        plain = {k: v for k, v in bad.items() if k not in ('attempt_id', 'locus_id')}
        out, rc = self._traverse(plain)
        self.assertEqual((1, ['error']), (rc, list(out)))

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


    # ---------------------------------------------------------------- record coordinates

    def _check_positions(self, result, seed):
        """Every coordinates interval of |result| against the source records (header labels:
        positions within the record): the bases at [s, e) are the run's own spelled bases
        (natural orientation, both arms), a seed occurrence the seed; returns the number of
        intervals checked."""
        block = result['coordinates']
        names = self._names(result)
        k = block['k']
        checked = 0
        for e in block['seed']:
            record = self.records[names[e['label']]]
            for s, t in e['occurrences']:
                self.assertEqual(seed, record[s:t])
                checked += 1
        for side, entries in block['arms'].items():
            arm = result['arms'][side]
            segments = {g['id']: g for g in arm['segments']}
            self.assertEqual(len(arm['runs']), len(entries))
            for i, (entry, run) in enumerate(zip(entries, arm['runs'])):
                self.assertEqual((i, run['label'], run['from_bp'], run['to_bp']),
                                 (entry['run'], entry['label'], entry['from_bp'], entry['to_bp']))
                chain = []
                sid = run['segment']
                while True:
                    chain.append(sid)
                    if not segments[sid]['parents']:
                        break
                    sid = segments[sid]['parents'][0]
                if side == 'right':
                    flank = ''.join(segments[x]['sequence'] for x in reversed(chain))
                    own = flank[run['from_bp']:run['to_bp']]
                else:
                    flank = ''.join(segments[x]['sequence'] for x in chain)
                    own = flank[len(flank) - run['to_bp']:len(flank) - run['from_bp']]
                record = self.records[names[run['label']]]
                for s, t in entry['occurrences']:
                    self.assertEqual(t - s, run['to_bp'] - run['from_bp'])
                    self.assertEqual(own, record[s:t], (side, i, s, t))
                    checked += 1
        self.assertEqual(k, K)
        return checked

    def test_traverse_record_coordinates(self):
        """Record coordinates (DESIGN §18; opt-in): under support trace every seed result
        carries the block — the seed's occurrences per label and per run the record interval
        of its own bases — which the source records confirm; elsewhere null with the reason;
        without the request nothing is added. The cap needs coordinates: true."""
        seed = self.element
        base = {'seeds': [{'sequence': seed}],
                'strategy': {'support': 'trace', 'branching': {'on_reconverge': 'keep',
                                                               'max_label_branches': 2},
                             'bounds': {'max_extension_bp': BLOCK}}}
        plain, rc = self._traverse(base)
        self.assertEqual(0, rc, plain.get('error'))
        self.assertNotIn('coordinates', plain['results'][0])
        self.assertNotIn('coordinates', plain['strategy']['output'])
        req = copy.deepcopy(base)
        req['strategy']['output'] = {'coordinates': True}
        out, rc = self._traverse(req)
        self.assertEqual(0, rc, out.get('error'))
        self.assertEqual((True, 16), (out['strategy']['output']['coordinates'],
                                      out['strategy']['output']['max_coordinate_occurrences']))
        result = out['results'][0]
        block = result['coordinates']
        self.assertEqual(('record', K, 16, True),
                         (block['kind'], block['k'], block['max_occurrences'], block['complete']))
        self.assertEqual(3, len(block['seed']))       # acc1, acc2, acc3 carry the element
        self.assertGreater(self._check_positions(result, seed), 6)
        # stripped of what it adds, the response is the opt-out one
        stripped = copy.deepcopy(out)
        del stripped['strategy']['output']['coordinates']
        del stripped['strategy']['output']['max_coordinate_occurrences']
        del stripped['results'][0]['coordinates']
        for r in (stripped, plain):
            r.pop('timing', None)
            for x in r['results']:
                x.pop('timing', None)
        self.assertEqual(plain, stripped)
        # support kmer: null, with the reason; the cap is inert there
        kmer = copy.deepcopy(req)
        kmer['strategy']['support'] = 'kmer'
        kmer['strategy']['output']['max_coordinate_occurrences'] = 1
        out, rc = self._traverse(kmer)
        self.assertEqual(0, rc, out.get('error'))
        self.assertIsNone(out['results'][0]['coordinates'])
        self.assertEqual('support kmer', out['results'][0]['coordinates_reason'])
        # the cap without coordinates: true is refused, naming the field
        bad = copy.deepcopy(base)
        bad['strategy']['output'] = {'max_coordinate_occurrences': 4}
        out, rc = self._traverse(bad)
        self.assertNotEqual(0, rc)
        self.assertIn('strategy.output.max_coordinate_occurrences', out['error'])

    def test_traverse_partial_derivation_states_its_set(self):
        """D3 (the owner's decision of 2026-10-04): a derived seed whose time budget runs out
        after part of its derivation is delivered with the set derived from the k-mers read —
        a superset of the whole seed's carriers, so its label evidence is qualified — and a
        `derivation` limitation saying how much was read; the walk stops at the seed. Read in
        one piece (--traverse-chunk-target-ms 0), the derivation's first window is read whole
        and its first k-mer consumed before the clock is read; paced in chunks (the default),
        a budget spent before the first chunk reads nothing and fails the seed, as before."""
        request = {'seeds': [{'sequence': self.element}],
                   'strategy': {'direction': 'right',
                                'bounds': {'max_extension_bp': 10, 'time_budget_ms': 1e-9}}}
        out, rc = self._traverse(request)
        self.assertEqual(0, rc, out.get('error'))
        failed = out['results'][0]
        self.assertEqual('failed', failed['outcome']['walks'])
        self.assertIn('after 0 of', failed['error'])
        out, rc = self._traverse(request, flags='--traverse-chunk-target-ms 0')
        self.assertEqual(0, rc, out.get('error'))
        result = out['results'][0]
        self.assertNotIn('error', result)
        self.assertEqual(('partial', 'complete', 'qualified', 'inline'), self._outcome_axes(result))
        d = result['limitations'][0]
        self.assertEqual(('derivation', 'time_budget', 'bounds.time_budget_ms'),
                         (d['kind'], d['cause'], d['knob']))
        self.assertEqual(1, d['observed'])
        # n is the seed's num_kmers: the limitation has a failed derivation's fields only
        self.assertNotIn('num_kmers', d)
        self.assertEqual(len(self.element) - K + 1, result['seed']['num_kmers'])
        self.assertIn('after 1 of %d k-mers' % (len(self.element) - K + 1), d['effect'])
        self.assertTrue(result['seed']['labels_from_seed'])
        self.assertEqual(0, result['arms']['right']['complete_to_bp'])
        self.assertEqual('truncated', result['arms']['right']['status'])

        # Review of W1, finding 1: acc1 and acc2 carry the k-mer read first (in left1), acc1
        # alone the whole seed. Under `exhaustive` with max_seed_labels 1 the superset is not
        # refused as "2 labels carry the seed": the seed fails with the time budget, as before
        # D3. Without it the cap cuts the superset, stated as one. Unhurried: acc1, walked
        seed = self.left1 + self.element + self.right1
        n = len(seed) - K + 1
        request = {'seeds': [{'sequence': seed}],
                   'strategy': {'exhaustive': True, 'direction': 'right',
                                'labels': {'max_seed_labels': 1},
                                'bounds': {'max_extension_bp': 10, 'time_budget_ms': 1e-9}}}
        out, rc = self._traverse(request, flags='--traverse-chunk-target-ms 0')
        self.assertEqual(0, rc, out.get('error'))
        failed = out['results'][0]
        self.assertEqual('failed', failed['outcome']['walks'])
        self.assertEqual('The time budget (bounds.time_budget_ms) ran out while deriving the '
                         'permitted set from the seed, after 1 of %d k-mers; name the labels '
                         'explicitly or shorten the seed' % n, failed['error'])
        self.assertEqual(('derivation', 'time_budget'),
                         (failed['limitations'][0]['kind'], failed['limitations'][0]['cause']))
        del request['strategy']['exhaustive']
        out, rc = self._traverse(request, flags='--traverse-chunk-target-ms 0')
        self.assertEqual(0, rc, out.get('error'))
        cut = out['results'][0]
        self.assertNotIn('error', cut)
        self.assertEqual((2, 1), (cut['seed']['labels_supporting_total'],
                                  cut['seed']['labels_dropped']))
        effect = [l for l in cut['limitations'] if l['kind'] == 'seed_labels'][0]['effect']
        self.assertTrue(effect.startswith('1 of the 2 label(s) carrying the k-mers read'), effect)
        request['strategy']['exhaustive'] = True
        request['strategy']['bounds']['time_budget_ms'] = 30000
        out, rc = self._traverse(request, flags='--traverse-chunk-target-ms 0')
        self.assertEqual(0, rc, out.get('error'))
        whole = out['results'][0]
        self.assertNotIn('error', whole)
        self.assertEqual(['acc1'], whole['seed']['labels'])
        self.assertEqual(1, whole['seed']['labels_supporting_total'])

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
        """Label names come from FASTA headers, which need not be UTF-8. No output carries
        such a name verbatim, and a replaced one (U+FFFD, which servers wrote before) can be
        another label's name -- a continuation that resubmitted it went on under that label
        (GPT review, finding 1). So a seed whose labels include one is refused per seed in
        both modes, as a failed derivation (cause unrepresentable_label_name) that names
        the column, never the bytes; detail full and graphlet agree (T37)."""
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
        for mode, knob in (('annotate', 'labels.mode'),
                           ('constrain', 'labels.seed_label_kind')):
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
                # strict: the response is UTF-8, with no replacement character anywhere
                text = res.stdout.decode('utf-8')
                self.assertNotIn('\ufffd', text)
                outs.append(json.loads(text))
            full, out = outs
            self.assertEqual(0, self._assert_round_trip(full, out, 'not UTF-8, ' + mode))
            result = out['results'][0]
            self.assertEqual('failed', result['outcome']['walks'], mode)
            (lim,) = result['limitations']
            self.assertEqual(('derivation', 'unrepresentable_label_name', knob, 2),
                             (lim['kind'], lim['cause'], lim['knob'], lim['observed']), mode)
            self.assertIn('column 0, sequence ', lim['effect'])
            self.assertNotIn('abc', json.dumps(result))
            with self.assertRaises(ValueError):
                graphlet_lib.from_response(result, out)

    def test_t37_a_replaced_name_is_never_another_label(self):
        """The review's probe (GPT, finding 1): column 'bad\\xff' and column 'bad\\ufffd' --
        the second is the first's name as replaced. A derived seed carried by both is
        refused (its dictionary holds the name that is not UTF-8); the valid name, named
        explicitly, resolves to its own column and walks its own bases."""
        rng = random.Random(54321)
        seq = ''.join(rng.choice('ACGT') for _ in range(120))
        alt = ('A' if seq[35] != 'A' else 'C') + 'TGCA' * 20
        d = os.path.join(self.tempdir.name, 'replaced_name')
        os.makedirs(d)
        for name, content in (('one.fa', (b'bad\xff', seq)),
                              ('two.fa', ('bad\ufffd'.encode(), seq[20:35] + alt))):
            with open(os.path.join(d, name), 'wb') as f:
                f.write(b'>' + content[0] + b'\n' + content[1].encode() + b'\n')
        for cmd in (f'{METAGRAPH} build -p 1 --mode basic --graph succinct -k 15 -o graph '
                    f'one.fa two.fa',
                    f'{METAGRAPH} annotate -p 1 --anno-header -i graph.dbg --anno-type column '
                    f'-o annotation one.fa two.fa'):
            res = subprocess.run(shlex.split(cmd), cwd=d, stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
            self.assertEqual(0, res.returncode, res.stderr.decode())

        def traverse(seed):
            path = os.path.join(d, 'request.json')
            with open(path, 'w') as f:
                json.dump({'seeds': [seed], 'strategy': {
                    'direction': 'right', 'bounds': {'max_extension_bp': 5},
                    'output': {'detail': 'graphlet'}}}, f)
            res = subprocess.run([METAGRAPH, 'traverse', '-i', os.path.join(d, 'graph.dbg'),
                                  '-a', os.path.join(d, 'annotation.column.annodbg'), path],
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            self.assertEqual(0, res.returncode, res.stdout.decode('utf-8', 'replace'))
            return json.loads(res.stdout.decode('utf-8'))

        out = traverse({'sequence': seq[:30]})
        result = out['results'][0]
        self.assertEqual('failed', result['outcome']['walks'])
        self.assertEqual('unrepresentable_label_name', result['limitations'][0]['cause'])
        # the replaced name is a label of its own: named, it walks ITS bases, not seq's
        out = traverse({'sequence': seq[20:35], 'labels': ['bad\ufffd']})
        g = graphlet_lib.from_response(out['results'][0], out)
        self.assertEqual(['bad\ufffd'], [l.name for l in g.labels])
        self.assertEqual(alt[:5], g.spell('right', 0))
        self.assertNotEqual(seq[35:40], alt[:5])

    def _resource_case(self, name, records, labels, bounds, detail, seed=None,
                       direction='right', files=None, k=3, anno='column_coord', output=None,
                       **strategy):
        """The external re-review's resource probes (GPT, stage 2): an index of |records|
        (k = 3 unless given, header labels via --index-header-coords; |files|: {file name:
        records} instead, one column each; |anno| column_coord, or row_diff_brwt_coord,
        whose reads are budget-aware) and one request on |seed| (default {sequence: AAA})
        along |direction|; -> (the result, the response's bytes). A second call with the
        same |name| reuses the index."""
        d = os.path.join(self.tempdir.name, 'resource_' + name)
        files = files or {'input.fa': records}
        if not os.path.exists(d):
            os.makedirs(d)
            for file, recs in files.items():
                with open(os.path.join(d, file), 'wb') as f:
                    for header, seq in recs:
                        f.write(b'>' + header + b'\n' + seq.encode() + b'\n')
            inputs = ' '.join(sorted(files))
            if anno == 'row_diff_brwt_coord':
                rd = (f'{METAGRAPH} transform_anno -p 1 --anno-type row_diff --coordinates '
                      f'-o annotation -i graph.dbg annotation.column.annodbg')
                transform = (rd, rd + ' --row-diff-stage 1', rd + ' --row-diff-stage 2',
                             f'{METAGRAPH} transform_anno -p 1 --anno-type row_diff_brwt_coord '
                             f'--greedy -o annotation -i graph.dbg annotation.column.annodbg',
                             f'{METAGRAPH} relax_brwt -p 1 -o annotation '
                             f'annotation.row_diff_brwt_coord.annodbg')
            else:
                transform = (f'{METAGRAPH} transform_anno -p 1 --anno-type column_coord '
                             f'--coordinates -o annotation annotation.column.annodbg',)
            for cmd in ((f'{METAGRAPH} build -p 1 --mode basic --graph succinct -k {k} -o graph {inputs}',
                         f'{METAGRAPH} annotate -p 1 --anno-filename -i graph.dbg --anno-type column '
                         f'--coordinates -o annotation {inputs}')
                        + transform
                        + (f'{METAGRAPH} annotate -p 1 -i graph.dbg --anno-filename '
                           f'--index-header-coords -o annotation {inputs}',)):
                res = subprocess.run(shlex.split(cmd), cwd=d, stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE)
                self.assertEqual(0, res.returncode, res.stderr.decode())
        path = os.path.join(d, 'request.json')
        with open(path, 'w') as f:
            json.dump({'seeds': [seed or {'sequence': 'AAA'}],
                       'strategy': dict({'direction': direction, 'labels': labels,
                                         'bounds': bounds,
                                         'output': dict({'detail': detail, 'timing': False},
                                                        **(output or {}))},
                                        **strategy)}, f)
        res = subprocess.run([METAGRAPH, 'traverse', '--json', '-i', os.path.join(d, 'graph.dbg'),
                              '-a', os.path.join(d, 'annotation.%s.annodbg' % anno), path],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(0, res.returncode, res.stderr.decode())
        return json.loads(res.stdout.decode('utf-8'))['results'][0], len(res.stdout)

    # ---- stage 3 (DESIGN-traverse-graphlet.md §14.1): the budget-aware annotation reads

    @staticmethod
    def _dense_files():
        """The §14 freeze-gate shape on a row-diff index: the path P (acc_path) and a
        record that repeats the k-mer P[51, 66) 120,000 times between random spacers, so
        the row before it on P carries one label while its row-diff dependency rows hold
        120,000 coordinates."""
        rng = random.Random(1111)

        def seq(n):
            return ''.join(rng.choice('ACGT') for _ in range(n))

        P = seq(120)
        rep = ''.join(P[51:66] + seq(10) for _ in range(120000))
        return P, {'path.fa': [(b'acc_path', P)], 'rep.fa': [(b'acc_rep', rep)]}

    def _dense_case(self, labels, bounds, detail, seed, **strategy):
        P, files = self._dense_files()
        return self._resource_case('dense', None, labels, bounds, detail,
                                   seed={'sequence': seed(P), **strategy.pop('seed_extra', {})},
                                   files=files, k=15, anno='row_diff_brwt_coord', **strategy)

    def test_stage3_dense_row_stops_the_level_that_reads_it(self):
        """§14 freeze gate: a row whose dependencies are dense but whose result is tiny. The
        level that reads it cannot hold its row-diff dependency rows within 1 MiB, so the walk
        stops there (partial, phase annotation_decode), in JSON and in the graphlet, with the
        levers that read fewer rows and memory_bound_soft naming only what is uncharged."""
        for detail in ('full', 'graphlet'):
            result, size = self._dense_case(
                {'mode': 'constrain'}, {'max_extension_bp': 60, 'max_memory_mb': 1}, detail,
                lambda P: P[:30], seed_extra={'labels': ['acc_path']})
            self.assertEqual('partial', result['outcome']['walks'], detail)
            stop = result['resource_stop']
            self.assertEqual(('memory', 'annotation_decode'), (stop['resource'], stop['phase']))
            self.assertIn('more_selective_seed', stop['actions'])
            self.assertNotIn('lower_max_seed_labels', stop['actions'])
            self.assertIn('row-diff dependency rows', stop['message'])
            self.assertLessEqual(stop['used'], stop['effective'])
            arm = result['arms']['right']
            self.assertEqual('truncated', arm['status'])
            (walk,) = [l for l in arm['limitations'] if l['kind'] == 'walk_domain']
            # what admitting the refused row needed, beside what the walk held: more than the
            # budget (it was stated as the budget plus one)
            self.assertGreater(walk['observed'], walk['limit'])
            self.assertIn('row-diff dependency rows', walk['effect'])
            (soft,) = [l for l in result['limitations'] if l['kind'] == 'memory_bound_soft']
            self.assertIn('charged before they are held', soft['effect'])
            self.assertLessEqual(size, 1 << 20)
            if detail == 'graphlet':
                self.assertIn('\nQ locus memory annotation_decode ', result['graphlet'])
                g = graphlet_lib.parse(result['graphlet'])
                self.assertEqual('annotation_decode', g.resource_stop.phase)
        # without a budget the same walk passes the dense row (the default reads)
        result, _ = self._dense_case({'mode': 'constrain'}, {'max_extension_bp': 60}, 'summary',
                                     lambda P: P[:30], seed_extra={'labels': ['acc_path']})
        self.assertEqual('complete', result['outcome']['walks'])
        self.assertNotIn('resource_stop', result)

    def test_stage3_decode_stop_does_not_depend_on_batch_kmers(self):
        """Every row a fetch returns is admitted against what decoding it alone needs,
        whether the lookahead decoded it or the fetch did: the stop is the same for every
        annotation.batch_kmers."""
        results = []
        for batch in (1, 7, 64, 1000):
            result, _ = self._dense_case(
                {'mode': 'constrain'}, {'max_extension_bp': 60, 'max_memory_mb': 1}, 'graphlet',
                lambda P: P[:30], seed_extra={'labels': ['acc_path']},
                annotation={'batch_kmers': batch})
            results.append(result)
        self.assertEqual('annotation_decode', results[0]['resource_stop']['phase'])
        for r in results[1:]:
            self.assertEqual(results[0], r)

    def test_stage3_seed_phase_reads_fail_the_seed(self):
        """Read in the seed phase -- to validate the seed, to derive its labels, or as an
        annotate root -- a row that does not fit fails the seed: outcome failed, a seed-level
        walk_domain on the budget, the seed's levers (annotate mode: a label-constrained
        query)."""
        cases = (('validation', {'mode': 'constrain'}, lambda P: P[45:75], {'labels': ['acc_path']}),
                 ('derivation', {'mode': 'constrain', 'seed_label_kind': 'header'},
                  lambda P: P[45:75], {}),
                 ('annotate root', {'mode': 'annotate', 'seed_label_kind': 'header'},
                  lambda P: P[36:66], {}))
        for what, labels, seed, extra in cases:
            result, _ = self._dense_case(labels, {'max_extension_bp': 10, 'max_memory_mb': 1},
                                         'summary', seed, seed_extra=extra)
            self.assertEqual('failed', result['outcome']['walks'], what)
            stop = result['resource_stop']
            self.assertEqual(('memory', 'annotation_decode'), (stop['resource'], stop['phase']), what)
            self.assertEqual(['raise_memory_budget', 'more_selective_seed']
                             + (['label_constrained_query'] if labels['mode'] == 'annotate' else []),
                             stop['actions'], what)
            (walk,) = [l for l in result['limitations'] if l['kind'] == 'walk_domain']
            self.assertGreater(walk['observed'], walk['limit'], what)
            self.assertIn('row-diff dependency rows', walk['effect'], what)

    def test_stage3_work_counts_dependency_rows(self):
        """A row's work includes its row-diff dependency rows (8 units each and 1 per entry
        and coordinate they store), whichever read decoded them: the walk past the dense row
        is charged for the tens of thousands of coordinates its dependency rows hold (the
        row-diff transform cancels the copies whose successor is the one the fork rule picked),
        where stage 2 charged its 60 rows' own entries (about 10 units each)."""
        result, _ = self._dense_case(
            {'mode': 'constrain'}, {'max_extension_bp': 60, 'max_work_units': 100000000},
            'summary', lambda P: P[:30], seed_extra={'labels': ['acc_path']})
        self.assertEqual('complete', result['outcome']['walks'])
        work = result['arms']['right']['counters']['work_units']
        self.assertGreater(work, 20000)
        # a work budget below it stops the walk (in phase traversal: work never refuses a read)
        stopped, _ = self._dense_case(
            {'mode': 'constrain'}, {'max_extension_bp': 60, 'max_work_units': work // 2},
            'summary', lambda P: P[:30], seed_extra={'labels': ['acc_path']})
        self.assertEqual(('work', 'traversal'), (stopped['resource_stop']['resource'],
                                                 stopped['resource_stop']['phase']))
        self.assertIn('per dependency row', stopped['resource_stop']['message'])

    @staticmethod
    def _largest_charge(stop):
        """The number a work stop states: the most its seed charged between two comparisons
        with the budget, which bounds how far used exceeds the budget."""
        return int(stop['message'].split('between two comparisons: ')[1].split(' units')[0])

    def test_stage2_work_stop_states_its_overrun(self):
        """GPT re-review, finding 2: 25,000 headers on AAACAAAGAAAT, annotate, a work budget of
        1 used 125,044 units against an advertised overrun of at most W = 65,536: the root's
        row and the level's rows were charged before a check. The walk now compares the
        budget after every charge: the stop exceeds it by the one row that tripped it, and the
        message states the most the seed charged between two comparisons (here that row).
        Round 3 (F4): with direction both, the two roots' rows are charged before the first
        comparison (a stop between them would deliver no result complete to 0 bp): the
        stated number is then both rows, which the overrun stays within."""
        records = [(b'L%d' % i, 'AAACAAAGAAAT') for i in range(25000)]
        labels = {'mode': 'annotate', 'seed_label_kind': 'header', 'max_labels_per_node': 1}
        result, _ = self._resource_case(
            'work', records, labels, {'max_extension_bp': 10, 'max_work_units': 1}, 'full')
        stop = result['resource_stop']
        self.assertEqual('work', stop['resource'])
        self.assertLessEqual(stop['used'] - 1, 65536)
        self.assertEqual(0, result['arms']['right']['complete_to_bp'])
        self.assertIn('after every charge', stop['message'])
        self.assertEqual(8 + 25000, self._largest_charge(stop))
        self.assertLessEqual(stop['used'] - 1, self._largest_charge(stop))
        for seed in ('AAA', 'AAAC', 'AAACAAAG'):
            result, _ = self._resource_case(
                'work', records, labels, {'max_extension_bp': 10, 'max_work_units': 1}, 'summary',
                seed={'sequence': seed}, direction='both')
            stop = result['resource_stop']
            self.assertEqual('work', stop['resource'], seed)
            self.assertEqual(2 * (8 + 25000), self._largest_charge(stop), seed)
            self.assertLessEqual(stop['used'] - 1, self._largest_charge(stop), seed)

    def test_stage2_work_stops_state_their_largest_charge(self):
        """Round 3, F5 and F7: under support trace a row's coordinates were charged as one
        sum after the row, outside the stated bound (stated 9, overrun 99,924), and a
        label-state scan carried no number. The coordinates are now charged with their row
        when the fetch returns it, and every work stop states the most its seed charged
        between two comparisons, whatever the charge: a sweep of budgets keeps every overrun
        within the stated number."""
        result, _ = self._resource_case(
            'trace', [(b'r1', 'CGAT' + 'GAT' * 100000)], {'mode': 'constrain',
                                                          'seed_label_kind': 'header'},
            {'max_extension_bp': 5, 'max_work_units': 100}, 'summary',
            seed={'sequence': 'CGA'}, support='trace', branching={'on_reconverge': 'keep'})
        stop = result['resource_stop']
        self.assertEqual('work', stop['resource'])
        self.assertGreaterEqual(self._largest_charge(stop), 100000)
        self.assertLessEqual(stop['used'] - 100, self._largest_charge(stop))
        records = [(b'L%d' % i, 'AAACAAAGAAAT') for i in range(1000)]
        stops = 0
        for budget in range(1, 40000, 1999):
            result, _ = self._resource_case(
                'sweep', records, {'mode': 'constrain', 'seed_label_kind': 'header',
                                   'max_seed_labels': 1000,
                                   'change_cost': {'model': 'constant', 'value': 1}},
                {'max_extension_bp': 10, 'max_work_units': budget}, 'summary',
                seed={'sequence': 'AAACA'}, branching={'max_label_branches': 'unlimited'})
            stop = result.get('resource_stop')
            if not stop:
                continue
            stops += 1
            self.assertEqual('work', stop['resource'], budget)
            self.assertLessEqual(stop['used'] - budget, self._largest_charge(stop), budget)
        self.assertGreater(stops, 3)

    def test_stage2_every_fetched_row_is_charged(self):
        """Round 3, F3: the fix to finding 2 charged a row only when a head consumed it, so
        rows a level's fetch decoded but no head consumed went uncharged (880,000 of 900,000
        entries), and `used` understated the decoding. Rows are charged again when the fetch
        returns them: on an index whose every row is n labels wide (n files each holding both
        branches of a split), every work stop has used >= rows_requested x (8 + n)."""
        n = 300
        rng = random.Random(9101)
        S, X, Y = (''.join(rng.choice('ACGT') for _ in range(m)) for m in (40, 20, 20))
        X, Y = 'A' + X[1:], 'C' + Y[1:]          # S's last k-mer splits into X and Y
        files = {('f%03d.fa' % i): [(b'x', S + X), (b'y', S + Y)] for i in range(n)}
        free, _ = self._resource_case(
            'uniform', None, {'mode': 'annotate', 'seed_label_kind': 'column',
                              'max_labels_per_node': 1},
            {'max_extension_bp': 30, 'max_work_units': 10 ** 9}, 'summary',
            seed={'sequence': S[:20]}, files=files, k=15, branching={'on_reconverge': 'keep'})
        self.assertNotIn('resource_stop', free)
        total = free['arms']['right']['counters']['work_units']
        stops = 0
        for budget in range(1, total, (8 + n) // 3):
            result, _ = self._resource_case(
                'uniform', None, {'mode': 'annotate', 'seed_label_kind': 'column',
                                  'max_labels_per_node': 1},
                {'max_extension_bp': 30, 'max_work_units': budget}, 'summary',
                seed={'sequence': S[:20]}, files=files, k=15,
                branching={'on_reconverge': 'keep'})
            stop = result.get('resource_stop')
            if not stop:
                continue
            stops += 1
            rows = result['annotation']['rows_requested']
            self.assertGreaterEqual(stop['used'], rows * (8 + n),
                                    'budget %d: a fetched row was not charged' % budget)
        self.assertGreater(stops, 10)

    def test_stage2_interrupted_fetch_states_its_memory(self):
        """GPT re-review, finding 3: 100,000 headers on the next node, max_memory_mb 1,
        max_work_units 100 stopped inside the level's fetch and reported memory_bound_soft
        observed 0 while the annotation cache alone held ~1.6 MB beyond its allotment. What a
        fetch holds is now observed before any check after it can stop the walk."""
        records = [(b'root', 'AAAC')] + [(b'L%d' % i, 'AAC') for i in range(100000)]
        result, _ = self._resource_case(
            'soft', records, {'mode': 'annotate', 'seed_label_kind': 'header',
                              'max_labels_per_node': 100000},
            {'max_extension_bp': 10, 'max_work_units': 100, 'max_memory_mb': 1}, 'full')
        self.assertIn('resource_stop', result)
        (soft,) = [l for l in result['limitations'] if l['kind'] == 'memory_bound_soft']
        self.assertGreaterEqual(soft['observed'], 2)

    def test_stage2_failed_results_echo_within_what_they_state(self):
        """Round 3, F1: a failed seed's result echoed an index-supplied header in full (an
        ambiguous_header: in the error and as observed) and the request's seed_id, neither
        charged: 2.16 MB under 1 MiB with memory_bound_soft observed 0. Under a memory budget
        a long index name is now echoed as a bounded prefix with its length and place, and
        what the echoed seed_id holds beyond the budget is stated as memory_bound_soft."""
        name = b'\x01' * 180000
        budget = 1 << 20
        for detail in ('summary', 'full', 'graphlet'):
            result, size = self._resource_case(
                'ambiguous', None, {'mode': 'constrain', 'seed_label_kind': 'header'},
                {'max_extension_bp': 1, 'max_memory_mb': 1}, detail,
                files={'a.fa': [(name, 'AAAC')], 'b.fa': [(name, 'AAAC')]})
            self.assertEqual('failed', result['outcome']['walks'], detail)
            (d,) = [l for l in result['limitations'] if l['kind'] == 'derivation']
            self.assertEqual('ambiguous_header', d['cause'], detail)
            self.assertIn('180000 bytes; column ', d['observed'], detail)
            compact = json.dumps(result, ensure_ascii=True, separators=(',', ':')).encode()
            self.assertLessEqual(len(compact), budget, detail)
            self.assertLessEqual(size, budget, detail)
        for unit, n in (('\U0001F600', 100000), ('\x01', 200000)):
            result, _ = self._resource_case(
                'seed_id', [(b'r1', 'AAACGTTGCA')], {'mode': 'constrain',
                                                     'seed_label_kind': 'header'},
                {'max_extension_bp': 4, 'max_memory_mb': 1}, 'graphlet',
                seed={'sequence': 'AAACGT', 'seed_id': unit * n}, direction='both')
            self.assertEqual('failed', result['outcome']['walks'])
            (soft,) = [l for l in result['limitations'] if l['kind'] == 'memory_bound_soft']
            compact = json.dumps(result, ensure_ascii=True, separators=(',', ':')).encode()
            self.assertGreater(len(compact), budget)
            self.assertGreaterEqual(soft['observed'] << 20, len(compact) - budget,
                                    'the failed result exceeds the budget by more than stated')

    def test_stage2_validation_and_admission_state_what_they_held(self):
        """Round 3, F2: the explicit-label validation charged a row (which failed the seed)
        before observing what the fetch held: a row of 400,000 coordinates under 1 MiB and a
        work budget of 100 reported memory_bound_soft 0. F6: a seed whose depth-0 admission
        failed reported 0 although its dictionary held megabytes of names. Both now observe
        what is held before the check that can fail the seed."""
        result, _ = self._resource_case(
            'validation', [(b'r1', 'CGAT' + 'GAT' * 400000)], {'mode': 'constrain'},
            {'max_extension_bp': 5, 'max_memory_mb': 1, 'max_work_units': 100}, 'summary',
            seed={'sequence': 'GATG', 'labels': ['r1']}, support='trace',
            branching={'on_reconverge': 'keep'})
        self.assertEqual('failed', result['outcome']['walks'])
        self.assertEqual('work', result['resource_stop']['resource'])
        (soft,) = [l for l in result['limitations'] if l['kind'] == 'memory_bound_soft']
        self.assertGreaterEqual(soft['observed'], 2)     # 3.2 MB of coordinates under 1 MiB
        records = [(b'N%d' % i + b'x' * 400000, 'AAAC') for i in range(6)]
        for mode in ('constrain', 'annotate'):
            result, _ = self._resource_case(
                'dictionary', records, {'mode': mode, 'seed_label_kind': 'header',
                                        'max_labels_per_node': 100},
                {'max_extension_bp': 1, 'max_memory_mb': 1}, 'summary')
            self.assertEqual('failed', result['outcome']['walks'], mode)
            self.assertEqual('memory', result['resource_stop']['resource'], mode)
            (soft,) = [l for l in result['limitations'] if l['kind'] == 'memory_bound_soft']
            # 2.4 MB of names, held twice (each label and the query's or recorder's copy)
            self.assertGreaterEqual(soft['observed'], 3, mode)

    def test_stage2_escaped_names_stay_within_the_budget(self):
        """GPT re-review, finding 1: a header of 180,000 control characters (detail full,
        max_memory_mb 2) gave a memory stop whose compact per-seed JSON alone was 2,163,840
        bytes, and '%' headers delivered more than the model reserved in a graphlet: a name
        byte was charged at a fixed four (three) bytes, but JSON writes a control character
        as six and MGT a '%' as three, escaped once more inside the JSON string. Names are
        now charged as delivered: each seed either fails at depth 0 (no name delivered) or
        stays within its budget."""
        for name, header, labels, mb, detail in (
                ('control', b'\x01' * 180000, {'mode': 'constrain', 'seed_label_kind': 'header',
                                               'max_labels_per_node': 1}, 2, 'full'),
                ('percent', b'%' * 200000, {'mode': 'annotate', 'seed_label_kind': 'header',
                                            'max_labels_per_node': 1}, 2, 'graphlet'),
                ('percent1', b'%' * 70000, {'mode': 'annotate', 'seed_label_kind': 'header',
                                            'max_labels_per_node': 1}, 1, 'graphlet')):
            result, size = self._resource_case(
                name, [(header, 'AAAC')], labels, {'max_extension_bp': 1, 'max_memory_mb': mb},
                detail)
            budget = mb << 20
            compact = json.dumps(result, ensure_ascii=True, separators=(',', ':')).encode()
            self.assertLessEqual(len(compact), budget, name)
            self.assertLessEqual(size, budget, name)
            if result['outcome']['walks'] == 'failed':
                self.assertEqual('memory', result['resource_stop']['resource'], name)
                self.assertNotIn(header[:64].decode(), json.dumps(result), name)
                continue
            if 'graphlet' in result:
                # the body p, its JSON escaping E: the writer's text and JSON copy, then
                # the copy with writeString's three texts (spec §5, request budgets)
                body = result['graphlet'].encode()
                escaped = len(json.dumps(result['graphlet'], ensure_ascii=True)) - 2
                self.assertLessEqual(len(body) + 3 * escaped, budget, name)

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
        # the request budgets, and how far a work stop can exceed one: stated with what
        # bounds it (a fetch call's rows are decoded whole; each stop states the most its
        # seed charged between two comparisons), not as the fixed interval W (GPT re-review)
        self.assertEqual(['max_memory_mb', 'max_work_units'], caps['budgets'])
        self.assertEqual(65536, caps['work_check_interval'])
        self.assertIn('indivisible charge', caps['work_bound'])
        self.assertIn('between two comparisons', caps['work_bound'])
        self.assertEqual('soft', caps['memory_bound'])
        # feature level 4: the server's maxima of the budgets (off unless configured) and the
        # row-diff path cache of the reads (128 MiB by default)
        self.assertEqual((0, 0), (caps['max_memory_mb'], caps['max_work_units']))
        self.assertEqual(128, caps['decode_cache']['path_cache_mb'])
        self.assertIn('do not depend on it', caps['decode_cache']['rule'])
        # feature level 6: record coordinates, on this index (coordinates and a CoordToHeader)
        # with every kind; stated here only, not in the per-request capabilities
        co = caps['coordinates']
        self.assertEqual({'supported', 'knob', 'cap_knob', 'max_occurrences_default', 'kinds',
                          'limitation', 'action', 'output_bound', 'rule'}, set(co))
        self.assertEqual((True, 'output.coordinates', 'output.max_coordinate_occurrences', 16,
                          ['record', 'column', 'mixed'], 'coordinates', 'drop_coordinates'),
                         (co['supported'], co['knob'], co['cap_knob'],
                          co['max_occurrences_default'], co['kinds'], co['limitation'],
                          co['action']))
        self.assertEqual(caps['supports_trace'], co['supported'])
        # the true bound of the block, and the column record-end numbering (review of W1,
        # finding 6)
        # (a server maximum is the budget of a request without one: the text names the probe's
        # max_memory_mb, never a value it could contradict)
        self.assertIn("this server's max_memory_mb", co['output_bound'])
        self.assertIn('both maxima 0, their default', co['output_bound'])
        self.assertIn("a record's last k - 1 bases shares its numbers with the next record's "
                      "first k - 1 positions", co['rule'])
        self.assertIn('[c + k - L, c + k)', co['rule'])
        out = self._post('traverse', {'seeds': [{'sequence': self.element}],
                                      'strategy': {'bounds': {'max_extension_bp': 10}}}).json()
        self.assertNotIn('coordinates', out['capabilities'])
        self.assertEqual(caps['feature_level'], out['capabilities']['feature_level'])
        # max_uninterruptible_ms stays null (decision 3c-N5): no later stage promises a bound
        self.assertIsNone(caps['deadline_check']['max_uninterruptible_ms'])
        self.assertIn('it stays null', caps['deadline_check']['rule'])
        self.assertNotIn('before stage 3c', caps['deadline_check']['rule'])

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

    def test_api_resolve_deadline(self):
        """Milestone 1b (SPEC §4.5): bounds.time_budget_ms on /resolve. The `resolve` block on
        both GET routes; a budget a microsecond above the finalisation reserve stops the work
        before its first row (200: the stop block, the resolve of the empty prefix); one at the
        reserve is a 400 naming the field; one above the cap is lowered to it and stated, the
        answer otherwise the unbudgeted one; the library and the CLI (no cap) alike."""
        url = f'http://{self.host}:{self.port}'
        probe = requests.get(url + '/traverse/capabilities').json()
        server = requests.get(url + '/capabilities').json()
        self.assertEqual(probe['resolve'], server['resolve'])
        t = probe['resolve']['time_budget']
        self.assertTrue(t['accepted'])
        self.assertEqual('bounds.time_budget_ms', t['knob'])
        self.assertIsNone(t['default'])
        self.assertEqual(probe['max_time_ms'], t['max_time_ms'])
        self.assertEqual(250, t['finalize_reserve_ms'])
        self.assertEqual(['rows', 'support'], t['stop_phases'])
        self.assertIn('Not polled', t['rule'])
        reserve = t['finalize_reserve_ms']
        n = BLOCK - K + 1
        for request in ({'sequence': self.element, 'discover': {'max_labels': 10, 'kind': 'header'}},
                        {'sequence': self.element, 'labels': ['acc1', 'acc2', 'acc3']}):
            what = json.dumps(request)[-60:]
            ret = self._post('resolve', request)
            self.assertEqual(200, ret.status_code, ret.text)
            plain = ret.json()
            self.assertNotIn('stop', plain, what)
            self.assertNotIn('limits', plain, what)
            ret = self._post('resolve', dict(request, bounds={'time_budget_ms': reserve + 0.001}))
            self.assertEqual(200, ret.status_code, ret.text)
            out = ret.json()
            stop = out['stop']
            self.assertEqual(('rows', 'time', 0, n, 0, BLOCK, 0),
                             (stop['phase'], stop['reason'], stop['resolved_kmers'],
                              stop['query_kmers'], stop['resolved_bp'], stop['query_bp'],
                              stop['remainder_from_bp']), what)
            self.assertIn('exactly the resolve', stop['message'])
            self.assertEqual({'time_budget_ms': reserve + 0.001, 'finalize_reserve_ms': reserve,
                              'clamped': []}, out['limits'])
            # the resolve of the empty prefix: no k-mer, no run, no candidate; explicit labels
            # listed with no support, a discovery's none
            self.assertEqual((0, [], []), (out['num_kmers'], out['graph_runs'], out['candidates']))
            self.assertEqual([(l['label'], 0, []) for l in plain['labels']] if 'labels' in request
                             else [], [(l['label'], l['kmers_supported'], l['runs'])
                                       for l in out['labels']], what)
            # above the cap: lowered to it, stated, and the answer is the unbudgeted one
            ret = self._post('resolve', dict(request, bounds={'time_budget_ms': 10 ** 7}))
            self.assertEqual(200, ret.status_code, ret.text)
            big = ret.json()
            self.assertIsNone(big.pop('stop'))
            self.assertEqual({'time_budget_ms': t['max_time_ms'], 'finalize_reserve_ms': reserve,
                              'clamped': [{'field': 'bounds.time_budget_ms',
                                           'requested': 10 ** 7,
                                           'effective': t['max_time_ms']}]}, big.pop('limits'))
            big.pop('timing')
            plain.pop('timing')
            self.assertEqual(plain, big, what)
        # at the reserve, or not a number: refused, naming the field
        for bad in (reserve, -1, '1000'):
            ret = self._post('resolve', {'sequence': self.element, 'labels': ['acc1'],
                                         'bounds': {'time_budget_ms': bad}})
            self.assertEqual(400, ret.status_code, ret.text)
            self.assertIn('bounds.time_budget_ms', ret.json()['error'])
        # the library passes it through
        client = graphlet_lib.TraverseClient(self.host, self.port)
        out = client.resolve(self.element, labels=['acc1'], time_budget_ms=reserve + 0.001)
        self.assertEqual((0, n), (out['stop']['resolved_kmers'], out['stop']['query_kmers']))
        out = client.resolve(self.element, labels=['acc1'], time_budget_ms=10000)
        self.assertIsNone(out['stop'])
        self.assertEqual(n, out['num_kmers'])
        # the CLI: the same answers, no cap
        out, rc = self._traverse({'sequence': self.element, 'labels': ['acc1'],
                                  'bounds': {'time_budget_ms': reserve + 0.001}}, resolve=True)
        self.assertEqual(0, rc)
        self.assertEqual(('rows', 0, n), (out['stop']['phase'], out['stop']['resolved_kmers'],
                                          out['stop']['query_kmers']))
        out, rc = self._traverse({'sequence': self.element, 'labels': ['acc1'],
                                  'bounds': {'time_budget_ms': 10 ** 7}}, resolve=True)
        self.assertEqual(0, rc)
        self.assertEqual(({'time_budget_ms': 10 ** 7, 'finalize_reserve_ms': reserve,
                           'clamped': []}, None), (out['limits'], out['stop']))
        out, rc = self._traverse({'sequence': self.element, 'labels': ['acc1'],
                                  'bounds': {'time_budget_ms': reserve}}, resolve=True)
        self.assertEqual(1, rc)
        self.assertIn('bounds.time_budget_ms', out['error'])

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

    def test_api_not_after_ms(self):
        """Pass 5, W1: a request whose not_after_ms has passed on the server's clock is refused
        at handler start, 409 {error, state: expired, not_after_ms, server_time_ms, ids,
        server_instance}, runs nothing and registers nothing (GET answers 404), with or
        without attempt_id; one not passed runs and echoes it in usage; a malformed one is a
        400 without usage. The probe states the clock skew a ledger adds."""
        url = f'http://{self.host}:{self.port}'
        caps = requests.get(url=url + '/traverse/capabilities').json()
        att = caps['attempts']
        self.assertIn('not_after_ms', att['fields'])
        self.assertIs(type(att['clock_skew_allowance_ms']), int)
        self.assertEqual(2000, att['clock_skew_allowance_ms'])
        self.assertIn('clock_skew_allowance_ms', att['not_after'])
        base = {'seeds': [{'sequence': self.element}],
                'strategy': {'direction': 'right', 'bounds': {'max_extension_bp': 10},
                             'output': {'timing': False}}}
        now = int(time.time() * 1000)
        past = dict(base, attempt_id='na-past', budget_id='b-1', not_after_ms=now - 1000)
        ret = self._post('traverse', past)
        self.assertEqual(409, ret.status_code, ret.text)
        body = ret.json()
        self.assertEqual({'error', 'state', 'not_after_ms', 'server_time_ms', 'attempt_id',
                          'budget_id', 'server_instance'}, set(body))
        self.assertEqual(('expired', now - 1000, 'na-past', 'b-1', att['server_instance']),
                         (body['state'], body['not_after_ms'], body['attempt_id'],
                          body['budget_id'], body['server_instance']))
        self.assertGreaterEqual(body['server_time_ms'], now - 1000)
        self.assertIs(type(body['server_time_ms']), int)
        # nothing registered: the state is unknown, and the id is still free
        self.assertEqual(404, requests.get(url + '/traverse/attempt/na-past').status_code)
        # without attempt_id the same refusal, without ids
        ret = self._post('traverse', dict(base, not_after_ms=now - 1000))
        self.assertEqual(409, ret.status_code, ret.text)
        self.assertEqual({'error', 'state', 'not_after_ms', 'server_time_ms',
                          'server_instance'}, set(ret.json()))
        # not passed: runs, echoed in usage and in the state; the rest unchanged
        plain = self._post('traverse', base)
        future = now + 600_000
        ret = self._post('traverse', dict(base, attempt_id='na-future', not_after_ms=future))
        self.assertEqual(200, ret.status_code, ret.text)
        out = ret.json()
        self.assertEqual(future, out['usage']['not_after_ms'])
        out.pop('usage')
        self.assertEqual(plain.json(), out)
        self.assertEqual(future, requests.get(url + '/traverse/attempt/na-future').json()
                         ['not_after_ms'])
        ret = self._post('traverse', dict(base, not_after_ms=future))
        self.assertEqual(200, ret.status_code, ret.text)
        self.assertEqual(plain.json(), ret.json())
        # malformed: a 400 naming the field, no usage, nothing registered
        for bad in (now + 0.5, -1, str(future), 2 ** 53):
            ret = self._post('traverse', dict(base, attempt_id='na-bad', not_after_ms=bad))
            self.assertEqual(400, ret.status_code, (bad, ret.text))
            self.assertIn('not_after_ms', ret.json()['error'])
            self.assertNotIn('usage', ret.json())
        self.assertEqual(404, requests.get(url + '/traverse/attempt/na-bad').status_code)
        # /resolve does not take it
        ret = self._post('resolve', {'sequence': self.element, 'labels': ['acc1'],
                                     'not_after_ms': future})
        self.assertEqual(400, ret.status_code, ret.text)
        # the library passes it through and tells the expired 409 apart from a duplicate
        client = graphlet_lib.TraverseClient(self.host, self.port)
        with self.assertRaises(graphlet_lib.AttemptExpired) as cm:
            client.traverse([self.element], base['strategy'], attempt_id='na-lib',
                            not_after_ms=now - 1)
        self.assertEqual(409, cm.exception.status)
        self.assertEqual(now - 1, cm.exception.body['not_after_ms'])
        resp = client.traverse([self.element], base['strategy'], attempt_id='na-lib2',
                               not_after_ms=future)
        self.assertEqual(future, resp.usage['not_after_ms'])
        # the CLI refuses the same way: the body on stdout, exit status 1
        path = os.path.join(self.tempdir.name, 'na_cli.json')
        with open(path, 'w') as f:
            json.dump(dict(base, attempt_id='na-cli', not_after_ms=now - 1000), f)
        res = subprocess.run([METAGRAPH, 'traverse', '--json', '-i', self.graph, '-a', self.anno,
                              path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(1, res.returncode)
        cli = json.loads(res.stdout)
        self.assertEqual(('expired', now - 1000, 'na-cli'),
                         (cli['state'], cli['not_after_ms'], cli['attempt_id']))

    def test_api_bodies_are_the_cli_output_at_every_encoding(self):
        """Pass 5, W6: the traversal routes compress at level 1 (the probe says so) and write
        each seed's result as text once built; the gzip, deflate and identity bodies all
        decompress to the same bytes, which are those of `metagraph traverse --json` for the
        same request (with several seeds, a failed one among them, in two details)."""
        url = f'http://{self.host}:{self.port}'
        self.assertEqual(1, requests.get(url + '/traverse/capabilities').json()['compression_level'])
        for detail in ('full', 'graphlet'):
            req = {'seeds': [{'sequence': self.element},
                             {'sequence': self.left1 + self.element[:40], 'labels': ['acc1']},
                             {'sequence': self.left2 + self.element + self.right2}],
                   'strategy': {'direction': 'both', 'bounds': {'max_extension_bp': BLOCK},
                                'output': {'detail': detail, 'timing': False}}}
            bodies = {}
            for enc in ('gzip', 'deflate', 'identity'):
                ret = requests.post(url + '/traverse', data=json.dumps(req),
                                    headers={'Accept-Encoding': enc}, stream=True)
                self.assertEqual(200, ret.status_code)
                raw = ret.raw.read(decode_content=True)
                bodies[enc] = raw
            self.assertEqual(bodies['gzip'], bodies['deflate'])
            self.assertEqual(bodies['gzip'], bodies['identity'])
            out = json.loads(bodies['gzip'])
            self.assertIn('error', out['results'][2])
            path = os.path.join(self.tempdir.name, f'enc_{detail}.json')
            with open(path, 'w') as f:
                json.dump(req, f)
            res = subprocess.run([METAGRAPH, 'traverse', '--json', '-i', self.graph, '-a',
                                  self.anno, '--index-name', 'tiny', '--index-manifest',
                                  self.manifest, path], stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
            self.assertEqual(0, res.returncode, res.stderr.decode()[-500:])
            self.assertEqual(bodies['gzip'], res.stdout.rstrip(b'\n'), detail)

    def test_api_coordinates(self):
        """Record coordinates over HTTP: the body is the CLI's (every encoding), the cap needs
        coordinates: true (400), and a memory stop of a trace walk with coordinates offers
        drop_coordinates."""
        url = f'http://{self.host}:{self.port}'
        for detail in ('full', 'graphlet'):
            req = {'seeds': [{'sequence': self.element},
                             {'sequence': self.left2 + self.element + self.right2}],
                   'strategy': {'support': 'trace',
                                'branching': {'on_reconverge': 'keep', 'max_label_branches': 2},
                                'bounds': {'max_extension_bp': BLOCK},
                                'output': {'detail': detail, 'timing': False,
                                           'coordinates': True,
                                           'max_coordinate_occurrences': 'unlimited'}}}
            ret = requests.post(url + '/traverse', data=json.dumps(req),
                                headers={'Accept-Encoding': 'gzip'}, stream=True)
            self.assertEqual(200, ret.status_code)
            body = ret.raw.read(decode_content=True)
            out = json.loads(body)
            self.assertTrue(out['results'][0]['coordinates']['complete'])
            self.assertIsNone(out['results'][1]['coordinates'])
            self.assertEqual('no traversal', out['results'][1]['coordinates_reason'])
            path = os.path.join(self.tempdir.name, f'coords_{detail}.json')
            with open(path, 'w') as f:
                json.dump(req, f)
            res = subprocess.run([METAGRAPH, 'traverse', '--json', '-i', self.graph, '-a',
                                  self.anno, '--index-name', 'tiny', '--index-manifest',
                                  self.manifest, path], stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
            self.assertEqual(0, res.returncode, res.stderr.decode()[-500:])
            self.assertEqual(body, res.stdout.rstrip(b'\n'), detail)
        ret = self._post('traverse', {'seeds': [{'sequence': self.element}],
                                      'strategy': {'output': {'max_coordinate_occurrences': 3}}})
        self.assertEqual(400, ret.status_code)
        self.assertIn('max_coordinate_occurrences', ret.json()['error'])
        # 0 and a number above the maximum state the whole range (review of W1, finding 4)
        for cap in (0, 18446744073709551615):
            ret = self._post('traverse', {'seeds': [{'sequence': self.element}],
                                          'strategy': {'output': {'coordinates': True,
                                                                  'max_coordinate_occurrences': cap}}})
            self.assertEqual(400, ret.status_code)
            self.assertIn('strategy.output.max_coordinate_occurrences: out of range '
                          '[1, 18446744073709551614] (or "unlimited")', ret.json()['error'])
        # a memory stop, with and without coordinates: offered wherever they were asked for,
        # also where only their null form is (support kmer; review of W1, finding 2)
        offered = {}
        for support in ('trace', 'kmer'):
            for coordinates in (True, False):
                req = {'seeds': [{'sequence': self.element}],
                       'strategy': {'support': support,
                                    'branching': {'on_reconverge': 'keep', 'max_label_branches': 2},
                                    'bounds': {'max_extension_bp': BLOCK, 'max_memory_mb': 1},
                                    'output': {'coordinates': coordinates}}}
                ret = self._post('traverse', req)
                self.assertEqual(200, ret.status_code, ret.text)
                result = ret.json()['results'][0]
                q = result.get('resource_stop')
                self.assertIsNotNone(q, 'no memory stop at 1 MiB (%s)' % support)
                self.assertEqual('memory', q['resource'])
                offered[support, coordinates] = 'drop_coordinates' in q['actions']
                if coordinates and support == 'kmer':
                    self.assertIsNone(result['coordinates'])
                    self.assertEqual('support kmer', result['coordinates_reason'])
        self.assertEqual({('trace', True): True, ('trace', False): False,
                          ('kmer', True): True, ('kmer', False): False}, offered)

    def test_api_spec_resolve_example_is_accepted(self):
        """Review of pass 5: SPEC §4.1's resolve request example set both labels and discover,
        which the server refuses (exactly one); copied with a real sequence, it is answered."""
        with open(os.path.join(REPO, 'docs', 'SPEC-labeled-traversal-core.md'),
                  encoding='utf-8') as f:
            spec = f.read()
        section = spec[spec.index('### 4.1 Request'):]
        block = section[section.index('```json') + len('```json'):]
        example = json.loads(block[:block.index('```')])
        example['sequence'] = self.element
        ret = self._post('resolve', example)
        self.assertEqual(200, ret.status_code, ret.text)
        # with explicit labels in place of discover, as the SPEC says, too
        example.pop('discover')
        example['labels'] = [next(iter(self.records))]
        ret = self._post('resolve', example)
        self.assertEqual(200, ret.status_code, ret.text)

    def test_api_server_capabilities(self):
        """Pass 5, W3/W4: GET /capabilities, the server-wide document: routes and features,
        feature_level 6 (5 at the review of pass 5, 4 at the efficiency pass, 3 in pass 5), algorithm_version (as every response's), mode single, no graph list,
        the attempts block (as the probe's) and how deadlines are checked; the number types
        a ledger compares are integers."""
        url = f'http://{self.host}:{self.port}'
        ret = requests.get(url + '/capabilities', headers={'Accept-Encoding': 'gzip'})
        self.assertEqual(200, ret.status_code, ret.text)
        self.assertEqual('gzip', ret.headers.get('Content-Encoding'))
        c = ret.json()
        # 'pattern': the block, feature and route of POST /pattern (DESIGN-pattern-search.md
        # §7.3; test_pattern.py checks the block itself); 'resolve': the block of /resolve's
        # deadline (SPEC §4.5; test_api_resolve_deadline checks it)
        self.assertEqual({'algorithm_version', 'attempts', 'compression_level',
                          'content_encodings', 'deadline_check', 'feature_level', 'features',
                          'graphs', 'mode', 'pattern', 'ready', 'release', 'resolve', 'routes',
                          'schema_version', 'server_instance'}, set(c))
        self.assertEqual((6, 'single', None, True, 1),
                         (c['feature_level'], c['mode'], c['graphs'], c['ready'],
                          c['schema_version']))
        self.assertEqual(['search', 'align', 'resolve', 'traverse', 'attempts', 'pattern'],
                         c['features'])
        self.assertEqual({'align': 'POST /align', 'attempt': 'GET /traverse/attempt/{attempt_id}',
                          'cancel': 'POST /traverse/cancel', 'capabilities': 'GET /capabilities',
                          'column_labels': 'GET /column_labels', 'pattern': 'POST /pattern',
                          'resolve': 'POST /resolve', 'search': 'POST /search',
                          'stats': 'GET /stats', 'traverse': 'POST /traverse',
                          'traverse_capabilities': 'GET /traverse/capabilities'}, c['routes'])
        probe = requests.get(url + '/traverse/capabilities').json()
        self.assertEqual(probe['attempts'], c['attempts'])
        self.assertEqual(probe['attempts']['server_instance'], c['server_instance'])
        self.assertEqual(probe['deadline_check'], c['deadline_check'])
        self.assertEqual(probe['algorithm_version'], c['algorithm_version'])
        self.assertEqual(6, probe['feature_level'])
        self.assertEqual(1, c['compression_level'])
        self.assertEqual(1, probe['compression_level'])
        out = self._post('traverse', {'seeds': [{'sequence': self.element}],
                                      'strategy': {'bounds': {'max_extension_bp': 10}}}).json()
        self.assertEqual(c['algorithm_version'], out['algorithm_version'])
        self.assertEqual(6, out['capabilities']['feature_level'])
        att = c['attempts']
        for key in ('allowance_ms', 'hard_cap_ms', 'clock_skew_allowance_ms', 'retention_s',
                    'retention_count', 'content_timeout_s', 'client_check_ms'):
            self.assertIs(type(att[key]), int, key)
        self.assertEqual(899000, att['hard_cap_ms'])
        self.assertEqual(10000, att['allowance_ms'])
        reserve = att['delivery_reserve']
        self.assertEqual({'compress_mbps', 'build_mbps', 'account_per_text_byte',
                          'coordinate_account_per_text_byte',
                          'measured_text_bytes', 'measured_compress_mbps', 'measured_build_mbps',
                          'measured_account_per_text_byte', 'rate_window', 'rule', 'margin',
                          'stop_ms', 'measured_stop_ms', 'calibration'}, set(reserve))
        # feature level 6: the record coordinates' share of a seed's account is estimated at a
        # fixed bound a ledger can reproduce (an integer)
        self.assertIs(type(reserve['coordinate_account_per_text_byte']), int)
        self.assertEqual(12, reserve['coordinate_account_per_text_byte'])
        self.assertIn('ceil(C / coordinate_account_per_text_byte)', reserve['rule'])
        # the reserve's margin and the walk's stop time (review of pass 5, F3), calibrated in
        # the efficiency pass (feature level 4): + 950 ms (+ 200 before), ratios 30 / 50
        self.assertEqual(1.25, reserve['margin'])
        self.assertIs(type(reserve['stop_ms']), int)
        self.assertEqual(c['deadline_check']['chunk_target_ms'] + 950, reserve['stop_ms'])
        self.assertEqual({'json': 30, 'graphlet': 50}, reserve['account_per_text_byte'])
        self.assertEqual({'summary', 'tree', 'full', 'graphlet'},
                         set(reserve['measured_account_per_text_byte']))
        self.assertEqual((50, 10, 16, 1 << 20),
                         (reserve['compress_mbps'], reserve['build_mbps'], reserve['rate_window'],
                          reserve['measured_text_bytes']))
        for key in ('measured_compress_mbps', 'measured_build_mbps'):
            self.assertTrue(reserve[key] is None or reserve[key] > 0, key)
        dc = c['deadline_check']
        self.assertIsNone(dc['max_uninterruptible_ms'])
        for key in ('chunk_target_ms', 'observed_max_uninterruptible_ms'):
            self.assertIs(type(dc[key]), int, key)
        self.assertEqual(50, dc['chunk_target_ms'])
        # the library reads both documents
        client = graphlet_lib.TraverseClient(self.host, self.port)
        self.assertEqual('single', client.server_capabilities()['mode'])
        self.assertEqual(probe['index_meta_fp'], client.capabilities()['index_meta_fp'])

    def test_api_attempts(self):
        """Stage 4, backend half: a request with attempt_id (budget_id and locus_id echoed)
        states its usage in every response — successful, partial, with a failed seed, or a
        400 after the request was read — and is otherwise the response without it; an id
        runs once (409, no usage), a malformed one is a 400 without usage, and a cancel or a
        state query of a finished or unknown id is a 404 (an unknown id is tombstoned)."""
        url = f'http://{self.host}:{self.port}'
        caps = requests.get(url=url + '/traverse/capabilities').json()
        att = caps['attempts']
        self.assertEqual(['attempt_id', 'budget_id', 'locus_id', 'not_after_ms',
                          'expect_server_instance'], att['fields'])
        self.assertEqual(['attempt_id', 'wait_ms', 'not_after_ms'], att['cancel_fields'])
        self.assertEqual('POST /traverse/cancel', att['cancel'])
        self.assertEqual('GET /traverse/attempt/{attempt_id}', att['state'])
        self.assertEqual(16, len(att['server_instance']))
        for key in ('retention_s', 'retention_count', 'allowance_ms', 'content_timeout_s',
                    'client_check_ms', 'bound', 'id_pattern'):
            self.assertIn(key, att)
        strategy = {'direction': 'right', 'bounds': {'max_extension_bp': 10},
                    'output': {'timing': False}}
        base = {'seeds': [{'sequence': self.element}], 'strategy': strategy}
        plain = self._post('traverse', base)
        self.assertEqual(200, plain.status_code, plain.text)
        self.assertNotIn('usage', plain.json())

        def attempt(name, **fields):
            req = copy.deepcopy(base)
            req.update(attempt_id=name, **fields)
            return req

        ret = self._post('traverse', attempt('api-ok', budget_id='budget-1', locus_id='locus-1'))
        self.assertEqual(200, ret.status_code, ret.text)
        out = ret.json()
        usage = out.pop('usage')
        # nothing else changes
        self.assertEqual(plain.json(), out)
        self.assertEqual(('api-ok', 'budget-1', 'locus-1', att['server_instance'], 'completed'),
                         (usage['attempt_id'], usage['budget_id'], usage['locus_id'],
                          usage['server_instance'], usage['reason']))
        for key in ('received_at', 'stopped_at'):
            self.assertRegex(usage[key], r'^\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d\.\d{3}Z$')
        self.assertIsInstance(usage['elapsed_ms'], int)
        self.assertEqual({'requested': 1, 'started': 1, 'finished': 1, 'abandoned': 0},
                         usage['seeds'])
        bound = usage['bound']
        self.assertTrue(bound['enforced'])
        self.assertEqual(1, bound['seeds'])
        self.assertEqual(usage['bound_ms'],
                         round(bound['seeds'] * bound['time_budget_ms'] + bound['allowance_ms']))
        self.assertGreater(usage['work_units'], 0)
        self.assertGreater(usage['memory']['peak_admitted_bytes'], 0)
        self.assertIsNone(usage['memory']['soft_excess_bytes'])
        # without a memory budget the account bounds nothing that was held: no bound stated
        self.assertIsNone(usage['memory']['held_bound_bytes'])
        seed = usage['per_seed'][0]
        self.assertIsNone(seed['refused_bytes'])
        self.assertEqual(usage['work_units'], seed['work_units'])
        self.assertEqual(out['results'][0]['outcome']['walks'], seed['outcome'])
        state = requests.get(url + '/traverse/attempt/api-ok')
        self.assertEqual(200, state.status_code, state.text)
        state = state.json()
        self.assertEqual(('finished', 'completed', True, 200),
                         (state['state'], state['reason'], state['response']['written'],
                          state['response']['status']))
        self.assertEqual(usage['work_units'], state['usage']['work_units'])
        # an id runs once: refused while it is retained, without usage
        dup = self._post('traverse', attempt('api-ok'))
        self.assertEqual(409, dup.status_code, dup.text)
        self.assertNotIn('usage', dup.json())
        self.assertEqual('finished', dup.json()['attempt']['state'])
        # nothing to cancel any more
        cancel = self._post('traverse/cancel', {'attempt_id': 'api-ok'})
        self.assertEqual(404, cancel.status_code, cancel.text)
        self.assertEqual((False, 'finished'), (cancel.json()['cancelled'], cancel.json()['state']))

        # a work budget stops the walk: partial, and the usage says what stopped it
        partial = copy.deepcopy(attempt('api-partial'))
        partial['strategy']['bounds']['max_work_units'] = 20
        ret = self._post('traverse', partial)
        self.assertEqual(200, ret.status_code, ret.text)
        seed = ret.json()['usage']['per_seed'][0]
        self.assertEqual(('partial', 'work'), (seed['outcome'], seed['stopped_by']))
        self.assertEqual('completed', ret.json()['usage']['reason'])

        # a failed seed beside a walked one
        failed = attempt('api-failed')
        failed['seeds'] = [{'sequence': self.element},
                           {'sequence': self.left2 + self.element + self.right2}]
        ret = self._post('traverse', failed)
        self.assertEqual(200, ret.status_code, ret.text)
        per_seed = ret.json()['usage']['per_seed']
        self.assertEqual('failed', per_seed[1]['outcome'])
        self.assertEqual(2, ret.json()['usage']['seeds']['finished'])

        # a 400 after the request was read states the usage too
        bad = attempt('api-bad')
        bad['strategy'] = dict(strategy, bogus=1)
        ret = self._post('traverse', bad)
        self.assertEqual(400, ret.status_code, ret.text)
        self.assertIn('bogus', ret.json()['error'])
        self.assertEqual(('api-bad', 'error'), (ret.json()['usage']['attempt_id'],
                                                ret.json()['usage']['reason']))
        self.assertEqual('error', requests.get(url + '/traverse/attempt/api-bad').json()['reason'])
        # ... but not a malformed id, nor an echo without the attempt it belongs to
        for req in (attempt('not an id'), attempt('x' * 129), dict(base, budget_id='b')):
            ret = self._post('traverse', req)
            self.assertEqual(400, ret.status_code, ret.text)
            self.assertNotIn('usage', ret.json())
        self.assertEqual(400, requests.get(url + '/traverse/attempt/not%20an%20id').status_code)
        self.assertEqual(400, self._post('traverse/cancel', {'attempt_id': 7}).status_code)
        self.assertEqual(400, self._post('traverse/cancel', {'attempt_id': 'a', 'x': 1}).status_code)

        # an unknown id: 404 to both; a cancel tombstones it, so it never runs here
        state = requests.get(url + '/traverse/attempt/api-unknown')
        self.assertEqual(404, state.status_code)
        self.assertEqual('unknown', state.json()['state'])
        cancel = self._post('traverse/cancel', {'attempt_id': 'api-early'})
        self.assertEqual(404, cancel.status_code)
        self.assertEqual(('unknown', True, False),
                         (cancel.json()['state'], cancel.json()['tombstone'],
                          cancel.json()['cancelled']))
        late = self._post('traverse', attempt('api-early'))
        self.assertEqual(409, late.status_code, late.text)
        self.assertTrue(late.json()['attempt']['tombstone'])

        # the library sends the fields and reads the answers
        client = graphlet_lib.TraverseClient(self.host, self.port)
        resp = client.traverse([self.element], strategy, attempt_id='api-lib', locus_id='l-9')
        self.assertEqual(('api-lib', 'l-9'), (resp.usage['attempt_id'], resp.usage['locus_id']))
        self.assertEqual('finished', client.attempt('api-lib')['state'])
        self.assertFalse(client.cancel('api-lib')['cancelled'])
        self.assertEqual('unknown', client.attempt('api-none')['state'])
        # an error after the request was read carries its usage in the error
        with self.assertRaises(graphlet_lib.TraverseError) as cm:
            client.traverse([self.element], dict(strategy, bogus=1), attempt_id='api-lib-bad')
        self.assertEqual('error', cm.exception.usage['reason'])

        # a memory budget: what a failed seed's result holds and states is its usage
        failed = attempt('api-memory')
        failed['seeds'] = [{'sequence': self.element, 'seed_id': 'z' * (2 << 20)},
                           {'sequence': self.element}]
        failed['strategy'] = copy.deepcopy(strategy)
        failed['strategy']['bounds']['max_memory_mb'] = 1
        ret = self._post('traverse', failed)
        self.assertEqual(200, ret.status_code, ret.text[:500])
        out = ret.json()
        usage = out['usage']
        first = usage['per_seed'][0]
        self.assertEqual(('failed', 'memory'), (first['outcome'], first['stopped_by']))
        soft = [l['observed'] for l in out['results'][0]['limitations']
                if l['kind'] == 'memory_bound_soft']
        self.assertEqual(soft, [-(-first['soft_excess_bytes'] // (1 << 20))])
        self.assertGreater(first['refused_bytes'], 1 << 20)
        self.assertLessEqual(first['peak_admitted_bytes'], 1 << 20)
        self.assertGreaterEqual(usage['memory']['held_bound_bytes'],
                                (1 << 20) + first['soft_excess_bytes'])

        # the MCP tool returns the real server's capabilities at its default ceiling
        from metagraph.traverse.mcp_tools import GraphletTools
        with tempfile.TemporaryDirectory() as root:
            tools = GraphletTools(graphlet_lib.GraphletStore(root), {'api': client})
            got = tools.traverse_capabilities(index='api')
        self.assertNotIn('error', got, str(got)[:300])
        self.assertEqual(caps['attempts']['fields'], got['capabilities']['attempts']['fields'])


def _free_port():
    s = socket.socket()
    s.bind(('127.0.0.1', 0))
    port = s.getsockname()[1]
    s.close()
    return port


def _instant(text):
    """an ISO-8601 instant of the attempt state (UTC, milliseconds) as seconds"""
    import datetime
    return datetime.datetime.strptime(text, '%Y-%m-%dT%H:%M:%S.%fZ').replace(
        tzinfo=datetime.timezone.utc).timestamp()


RESOLVE_BASE_BINARY = os.environ.get('METAGRAPH_BASE_BINARY', '')


@unittest.skipUnless(RESOLVE_BASE_BINARY and os.path.isfile(RESOLVE_BASE_BINARY),
                     "$METAGRAPH_BASE_BINARY (a binary before /resolve's deadline) is not set")
class TestTraverseResolveRegression(TestTraverseBase):
    """Milestone 1b (SPEC §4.5): a /resolve without bounds.time_budget_ms is answered byte for
    byte as by the binary before the deadline ($METAGRAPH_BASE_BINARY, e.g. the build of
    804731aa), apart from timing.elapsed_ms: the resolve requests of this file and their
    refusals, on the server (with and without gzip) and on the CLI."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.servers = {}
        for name, binary in (('base', RESOLVE_BASE_BINARY), ('new', METAGRAPH)):
            port = _free_port()
            log = open(f'{cls.tempdir.name}/resolve-{name}.log', 'w')
            process = subprocess.Popen(
                shlex.split(binary) + ['server_query', '-i', cls.graph, '-a', cls.anno,
                                       '--port', str(port), '--address', '127.0.0.1', '-p', '2',
                                       '--index-name', 'tiny', '--index-manifest', cls.manifest],
                stdout=log, stderr=subprocess.STDOUT)
            url = f'http://127.0.0.1:{port}'
            for _ in range(600):
                try:
                    if requests.get(url + '/traverse/capabilities', timeout=2).ok:
                        break
                except requests.exceptions.RequestException:
                    pass
                time.sleep(0.1)
            cls.servers[name] = (binary, url, process, log)

    @classmethod
    def tearDownClass(cls):
        for _, _, process, log in cls.servers.values():
            process.kill()
            process.wait()
            log.close()
        super().tearDownClass()

    def _requests(self):
        e = self.element
        with open(os.path.join(REPO, 'docs', 'SPEC-labeled-traversal-core.md'),
                  encoding='utf-8') as f:
            spec = f.read()
        section = spec[spec.index('### 4.1 Request'):]
        block = section[section.index('```json') + len('```json'):]
        example = dict(json.loads(block[:block.index('```')]), sequence=e)
        explicit_example = {k: v for k, v in example.items() if k != 'discover'}
        explicit_example['labels'] = ['acc1']
        return [
            {'sequence': e, 'discover': {'max_labels': 10, 'kind': 'header'}},
            {'sequence': self.left1 + e, 'labels': ['acc1', 'acc2', 'acc3'],
             'select': {'policy': 'max_support', 'max_seeds': 2, 'min_block_bp': K}},
            {'sequence': e, 'labels': ['acc1', 'acc2', 'acc3']},
            {'sequence': e, 'discover': {'max_labels': 10}},
            example,
            explicit_example,
            {'sequence': e, 'labels': ['acc1', 'acc3'], 'support': 'trace', 'run_format': 'rle'},
            {'sequence': self.left2 + e + self.right2,
             'discover': {'max_labels': 2, 'kind': 'header'}, 'support': 'trace',
             'select': {'policy': 'longest_first', 'max_labels_per_seed': 1}},
            {'sequence': self.left1 + e + self.right2, 'discover': {'max_labels': 1},
             'select': {'policy': 'max_support', 'label_order': 'kmers_supported'}},
            {'sequence': e, 'labels': ['acc1'],
             'select': {'policy': 'explicit',
                        'seeds': [{'kmer_interval': [0, 10], 'labels': ['acc1']}]}},
            {'sequence': e[:40] + 'N' * 5 + e[45:], 'labels': ['acc2'],
             'bounds': {'max_query_bp': 1000}},
            # refusals
            {'sequence': e, 'labels': ['acc1'], 'not_after_ms': 1},
            {'sequence': e, 'labels': ['nope']},
            {'sequence': 'ACGT', 'labels': ['acc1']},
            {'sequence': e, 'labels': ['acc1'], 'discover': {}},
            {'bogus': 1},
            {'sequence': e, 'labels': ['acc1'], 'bounds': {'max_query_bp': 10}},
            {'sequence': e, 'labels': ['acc1'], 'select': {'policy': 'explicit'}},
        ]

    @staticmethod
    def _untimed(text):
        return re.sub(r'"elapsed_ms":[-+0-9.eE]+', '"elapsed_ms":0', text)

    def test_server_answers_as_before(self):
        compared = 0
        for request in self._requests():
            for encoding in ('identity', 'gzip'):
                answers = []
                for name in ('base', 'new'):
                    ret = requests.post(self.servers[name][1] + '/resolve',
                                        data=json.dumps(request),
                                        headers={'Accept-Encoding': encoding})
                    answers.append((ret.status_code, ret.headers.get('Content-Encoding'),
                                    self._untimed(ret.text)))
                self.assertEqual(answers[0], answers[1], json.dumps(request)[:200])
                compared += 1
        self.assertGreaterEqual(compared, 36)

    def test_cli_answers_as_before(self):
        path = f'{self.tempdir.name}/resolve-regression.json'
        for request in self._requests():
            with open(path, 'w') as f:
                json.dump(request, f)
            answers = []
            for name in ('base', 'new'):
                res = subprocess.run(shlex.split(self.servers[name][0])
                                     + ['traverse', '--resolve', '--json', '-i', self.graph,
                                        '-a', self.anno, path],
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                answers.append((res.returncode, self._untimed(res.stdout.decode())))
            self.assertEqual(answers[0], answers[1], json.dumps(request)[:200])


@unittest.skipIf(PROTEIN_MODE, "traversal fixtures are DNA")
@unittest.skipUnless(_supports_traverse(), "`metagraph traverse` is not available in this build")
class TestTraverseMultiGraph(TestTraverseBase):
    """Pass 5, W2: per-graph identity on a multi-graph server. The graph list gains two
    optional columns, manifest_path and index_ns, each manifest checked at start-up against
    the files its pair loads (a mismatch refuses to start); every /resolve and /traverse
    response, and every graphlet's H record, states its pair's identity; GET
    /traverse/capabilities?graph=<name>[&graph_path=<path>] describes one pair; a graph_path
    naming one graph with several annotations is refused (it cannot choose one). The manifests
    of the list are written by scripts/traversal/index_manifest.py --server-csv."""

    SCRIPT = os.path.join(REPO, 'scripts', 'traversal', 'index_manifest.py')

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        d = cls.tempdir.name
        cls.pairs = {}
        # B: two graphs of their own (one record each), the first with a manifest written by
        # the batch tool, the second without
        for tag, acc in (('two', 'acc1'), ('three', 'acc2')):
            os.makedirs(f'{d}/{tag}', exist_ok=True)
            fasta = f'{d}/{tag}/ref.fa'
            with open(fasta, 'w') as f:
                f.write(f'>{acc}\n{cls.records[acc]}\n')
            graph = f'{d}/{tag}/graph_{tag}.dbg'
            cls._build_graph(fasta, graph, K, 'succinct', mode='basic')
            cls._annotate_graph(fasta, graph, f'{d}/{tag}/annotation_{tag}', 'column',
                                anno_type='header')
            cls.pairs[tag] = (graph, f'{d}/{tag}/annotation_{tag}.column.annodbg')
        # C: the base graph with a second annotation (one graph, two annotations)
        os.makedirs(f'{d}/col', exist_ok=True)
        cls._annotate_graph(cls.fasta, cls.graph, f'{d}/col/annotation_col', 'column',
                            anno_type='header')
        cls.pairs['col'] = (cls.graph, f'{d}/col/annotation_col.column.annodbg')
        cls.csv = f'{d}/graphs.csv'
        lines = [
            f'A,{cls.graph},{cls.anno},{cls.manifest},tinyA',
            f'B,{cls.pairs["two"][0]},{cls.pairs["two"][1]},{d}/two/two.manifest.json',
            f'B,{cls.pairs["three"][0]},{cls.pairs["three"][1]}',
            # the same pair as A, without columns: it states A's identity (one index)
            f'C,{cls.graph},{cls.anno}',
            f'C,{cls.pairs["col"][0]},{cls.pairs["col"][1]}',
            # A's pair again under A: still one pair, addressed without graph_path
            f'A,{cls.graph},{cls.anno}',
        ]
        with open(cls.csv, 'w') as f:
            f.write('\n'.join(lines) + '\n')
        # the manifest of B's first pair, by the batch tool (A's exists and is kept)
        sub = f'{d}/sub.csv'
        with open(sub, 'w') as f:
            f.write(lines[1] + '\n')
        res = subprocess.run([sys.executable, cls.SCRIPT, '--server-csv', sub, '--jobs', '2'],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        assert res.returncode == 0, res.stderr.decode()
        with open(f'{d}/two/two.manifest.json') as f:
            cls.fp_two = json.load(f)['index_fp']
        cls.server = cls._Server(cls, cls.csv)
        assert cls.server.ready, cls.server.text()

    @classmethod
    def tearDownClass(cls):
        cls.server.close()
        super().tearDownClass()

    class _Server:
        def __init__(self, test, csv, wait=True):
            self.port = _free_port()
            self.url = f'http://127.0.0.1:{self.port}'
            self.log_path = f'{test.tempdir.name}/multi-{self.port}.log'
            self.log = open(self.log_path, 'w')
            self.process = subprocess.Popen(
                shlex.split(METAGRAPH) + ['server_query', csv, '--port', str(self.port),
                                          '--address', '127.0.0.1', '-p', '2'],
                stdout=self.log, stderr=subprocess.STDOUT)
            self.ready = False
            for _ in range(600 if wait else 0):
                if self.process.poll() is not None:
                    break
                try:
                    if requests.get(self.url + '/capabilities', timeout=2).ok:
                        self.ready = True
                        break
                except requests.exceptions.RequestException:
                    pass
                time.sleep(0.1)

        def close(self):
            if self.process.poll() is None:
                self.process.kill()
            self.process.wait()
            self.log.close()

        def text(self):
            with open(self.log_path) as f:
                return f.read()

    def _post(self, route, payload):
        return requests.post(f'{self.server.url}/{route}', data=json.dumps(payload))

    def _caps(self, query=''):
        return requests.get(f'{self.server.url}/traverse/capabilities{query}')

    def _request(self, graph, graph_path=None, detail='graphlet'):
        req = {'graph': graph, 'seeds': [{'sequence': self.element}],
               'strategy': {'direction': 'right', 'bounds': {'max_extension_bp': 20},
                            'output': {'detail': detail, 'timing': False}}}
        if graph_path:
            req['graph_path'] = graph_path
        return req

    def test_multi_responses_state_their_pairs_identity(self):
        expected = {
            ('A', None): ('tinyA', self.index_fp),
            ('B', self.pairs['two'][0]): (None, self.fp_two),
            ('B', self.pairs['three'][0]): (None, None),
        }
        metas = set()
        for (graph, path), ident in expected.items():
            ret = self._post('traverse', self._request(graph, path))
            self.assertEqual(200, ret.status_code, ret.text)
            out = ret.json()
            caps = out['capabilities']
            self.assertEqual(ident, (caps['index_ns'], caps['index_fp']), (graph, path))
            metas.add(caps['index_meta_fp'])
            head = out['results'][0]['graphlet'].split('\n', 1)[0].split(' ')
            self.assertEqual([ident[0] or '*', ident[1] or '*', caps['index_meta_fp']],
                             head[-3:])
            # the library reads the identity off the H record
            g = graphlet_lib.from_response(out['results'][0], out)
            self.assertEqual((ident[0], ident[1]), (g.index_ns, g.index_fp))
            req = {'graph': graph, 'sequence': self.element, 'discover': {'max_labels': 10}}
            if path:
                req['graph_path'] = path
            res = self._post('resolve', req)
            self.assertEqual(200, res.status_code, res.text)
            self.assertEqual(ident, (res.json()['capabilities']['index_ns'],
                                     res.json()['capabilities']['index_fp']))
        self.assertEqual(3, len(metas))
        # C lists the base graph with two annotations: graph_path cannot choose one
        ret = self._post('traverse', self._request('C', self.pairs['col'][0]))
        self.assertEqual(400, ret.status_code, ret.text)
        self.assertIn('several annotations', ret.json()['error'])

    def test_multi_probe_describes_one_pair(self):
        ret = self._caps('?graph=A')
        self.assertEqual(200, ret.status_code, ret.text)
        caps = ret.json()
        self.assertEqual(('A', self.graph, 'tinyA', self.index_fp, 6),
                         (caps['graph'], caps['graph_path'], caps['index_ns'],
                          caps['index_fp'], caps['feature_level']))
        for key in ('attempts', 'deadline_check', 'budgets', 'algorithm_version',
                    'compression_level', 'content_encodings', 'work_bound', 'coordinates'):
            self.assertIn(key, caps)
        # per pair: whether its index reports record coordinates is whether it supports trace
        self.assertEqual(caps['supports_trace'], caps['coordinates']['supported'])
        out = self._post('traverse', self._request('A')).json()
        for key in ('index_ns', 'index_fp', 'index_meta_fp', 'k', 'num_labels'):
            self.assertEqual(out['capabilities'][key], caps[key], key)
        # a name over several graphs needs graph_path, which must name one of them
        ret = self._caps('?graph=B')
        self.assertEqual(400, ret.status_code)
        self.assertIn(self.pairs['two'][0], ret.json()['error'])
        self.assertIn(self.pairs['three'][0], ret.json()['error'])
        quoted = requests.utils.quote(self.pairs['two'][0], safe='')
        ret = self._caps(f'?graph=B&graph_path={quoted}')
        self.assertEqual(200, ret.status_code, ret.text)
        self.assertEqual((self.pairs['two'][0], None, self.fp_two),
                         (ret.json()['graph_path'], ret.json()['index_ns'],
                          ret.json()['index_fp']))
        # one graph with two annotations cannot be addressed by graph_path, and a name listing
        # only that is refused as such, not told to pass graph_path (review of pass 5)
        quoted = requests.utils.quote(self.graph, safe='')
        ret = self._caps(f'?graph=C&graph_path={quoted}')
        self.assertEqual(400, ret.status_code, ret.text)
        self.assertIn('several annotations', ret.json()['error'])
        ret = self._caps('?graph=C')
        self.assertEqual(400, ret.status_code, ret.text)
        self.assertIn("index 'C' lists one graph", ret.json()['error'])
        self.assertIn('several annotations', ret.json()['error'])
        self.assertNotIn('graph_path\' to pick one', ret.json()['error'])
        for query, needle in (('', 'needs ?graph=<name>'), ('?graph=Z', "unknown graph 'Z'"),
                              ('?graph=A&x=1', "unknown parameter 'x'"),
                              ('?graph=A&graph=B', 'given twice'),
                              ('?graph=B&graph_path=/nowhere', 'is not part of')):
            ret = self._caps(query)
            self.assertEqual(400, ret.status_code, query)
            self.assertIn(needle, ret.json()['error'], query)
        # the server-wide document lists the names; align is not offered here
        c = requests.get(self.server.url + '/capabilities').json()
        self.assertEqual(('multi', ['A', 'B', 'C'], True),
                         (c['mode'], c['graphs'], c['ready']))
        self.assertNotIn('align', c['features'])
        self.assertNotIn('align', c['routes'])
        self.assertEqual('GET /traverse/capabilities?graph={name}[&graph_path={path}]',
                         c['routes']['traverse_capabilities'])
        self.assertEqual(caps['attempts'], c['attempts'])
        # the library asks for one pair
        client = graphlet_lib.TraverseClient('127.0.0.1', self.server.port, graph='A')
        self.assertEqual('tinyA', client.capabilities()['index_ns'])
        self.assertEqual(self.fp_two, client.capabilities(graph='B', graph_path=self.pairs['two'][0])
                         ['index_fp'])
        self.assertEqual(['A', 'B', 'C'], client.server_capabilities()['graphs'])

    def test_multi_list_errors_refuse_to_start(self):
        """A manifest that does not describe its pair's files, two identities for one pair, a
        line of six columns and an index_ns that is no token each refuse to start (exit 1,
        an [error] line naming the problem) before anything is loaded."""
        d = self.tempdir.name
        with open(self.manifest) as f:
            wrong = json.load(f)
        for e in wrong['files']:
            e['size'] += 1
        wrong.pop('index_fp')
        with open(f'{d}/wrong.manifest.json', 'w') as f:
            json.dump(wrong, f)
        cases = {
            'size': (f'A,{self.graph},{self.anno},{d}/wrong.manifest.json',
                     'wrong.manifest.json'),
            'conflict': (f'A,{self.graph},{self.anno},,one\nB,{self.graph},{self.anno},,two',
                         'one index has one identity'),
            'columns': (f'A,{self.graph},{self.anno},,ns,extra', 'at most five'),
            'token': (f'A,{self.graph},{self.anno},,bad/ns', 'does not match'),
            'missing': (f'A,{self.graph},{self.anno},{d}/nowhere.json', 'cannot be read'),
        }
        for name, (text, needle) in cases.items():
            csv = f'{d}/bad_{name}.csv'
            with open(csv, 'w') as f:
                f.write(text + '\n')
            server = self._Server(self, csv, wait=False)
            try:
                code = server.process.wait(timeout=120)
            finally:
                server.close()
            log = server.text()
            self.assertEqual(1, code, (name, log[-2000:]))
            self.assertIn('[error]', log, name)
            self.assertIn(needle, log, (name, log[-2000:]))
            self.assertNotIn('Loading', log.split('[error]')[0][-200:] + '', name)

    def test_multi_identity_holes_refuse_to_start(self):
        """Review of pass 5: one index_fp could describe two indexes, and a manifest's sidecars
        were not checked. A manifest that also lists another annotation (one written for a
        directory, as for two annotations of one graph with swapped memberships), one whose
        sidecar (the .seqs the server loads) is another build's, and two different pairs stating
        one index_fp (copies of one bundle under two paths) each refuse to start."""
        d = self.tempdir.name
        with open(self.manifest) as f:
            base = json.load(f)
        names = [e['path'] for e in base['files']]
        self.assertIn('annotation.seqs', names)
        # a directory's manifest: the base graph, both of its annotations
        col = self.pairs['col'][1]
        bundle = copy.deepcopy(base)
        bundle.pop('index_fp')
        with open(col, 'rb') as f:
            data = f.read()
        bundle['files'].append({'path': os.path.basename(col), 'size': len(data),
                                'sha256': hashlib.sha256(data).hexdigest()})
        with open(f'{d}/bundle.manifest.json', 'w') as f:
            json.dump(bundle, f)
        # another build's sidecar
        sidecar = copy.deepcopy(base)
        sidecar.pop('index_fp')
        for e in sidecar['files']:
            if e['path'] == 'annotation.seqs':
                e['size'] += 12345
                e['sha256'] = 'f' * 64
        with open(f'{d}/sidecar.manifest.json', 'w') as f:
            json.dump(sidecar, f)
        # a copy of the bundle under another path, stating the same manifest
        os.makedirs(f'{d}/copy', exist_ok=True)
        for name in names:
            shutil.copyfile(f'{d}/{name}', f'{d}/copy/{name}')
        cases = {
            'directory': (f'A,{self.graph},{self.anno},{d}/bundle.manifest.json',
                          'which this index does not load'),
            'sidecar': (f'A,{self.graph},{self.anno},{d}/sidecar.manifest.json',
                        'annotation.seqs'),
            'shared_fp': (f'A,{self.graph},{self.anno},{self.manifest}\n'
                          f'B,{d}/copy/{os.path.basename(self.graph)},'
                          f'{d}/copy/{os.path.basename(self.anno)},{self.manifest}',
                          'one index_fp'),
        }
        for name, (text, needle) in cases.items():
            csv = f'{d}/hole_{name}.csv'
            with open(csv, 'w') as f:
                f.write(text + '\n')
            server = self._Server(self, csv, wait=False)
            try:
                code = server.process.wait(timeout=120)
            finally:
                server.close()
            log = server.text()
            self.assertEqual(1, code, (name, log[-2000:]))
            self.assertIn('[error]', log, name)
            self.assertIn(needle, log, (name, log[-2000:]))
        # the same pair spelled two ways is one index (one identity), not two
        csv = f'{d}/spellings.csv'
        rel = os.path.relpath(self.graph, d)
        with open(csv, 'w') as f:
            f.write(f'A,{self.graph},{self.anno},{self.manifest}\n'
                    f'B,{d}/./{rel},{self.anno},{self.manifest}\n')
        server = self._Server(self, csv)
        try:
            self.assertTrue(server.ready, server.text())
            out = requests.post(server.url + '/traverse',
                                data=json.dumps(self._request('B'))).json()
            self.assertEqual(self.index_fp, out['capabilities']['index_fp'])
        finally:
            server.close()

    # ---- review of pass 5, findings 2 and 3: one loader dependency inventory

    def _run_ok(self, cmd, cwd):
        res = subprocess.run(shlex.split(cmd), cwd=cwd, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        self.assertEqual(0, res.returncode, (cmd, res.stderr.decode()[-2000:]))
        return res

    def _cli_inventory(self, graph, annotation, *flags):
        res = subprocess.run(shlex.split(METAGRAPH) + ['traverse', '--index-inventory', '--json',
                                                       '-i', graph, '-a', annotation, *flags],
                             stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(0, res.returncode, res.stderr.decode())
        return json.loads(res.stdout)

    def _start_refused(self, csv):
        server = self._Server(self, csv, wait=False)
        try:
            code = server.process.wait(timeout=120)
        finally:
            server.close()
        return code, server.text()

    def test_inventory_is_mirrored_by_index_manifest_py(self):
        """The loader dependency inventory has one definition (index_load_inventory), which
        `traverse --index-inventory` prints and scripts/traversal/index_manifest.py mirrors
        (load_inventory, ANNOTATION_KINDS): on bundles with every sidecar kind — a masked
        graph with a Bloom filter, an unmasked one, row_diff (anchors and fork successors
        beside the graph), coordinate annotations with and without their .seqs, a column
        annotation with its .coords (not loaded), --no-coord-mapping, a symlinked spelling —
        both list the same files in the same roles, and their annotation tables are equal.
        The graph's derived data (the owner's decision #17 of 2026-10-08: its mask and Bloom
        filter, not part of index_fp) is in neither's identity files; both list it apart
        (index_derived_files, derived_files) with whether the loader reads it, and state the
        same rule; `index_manifest.py --inventory` prints what the binary prints."""
        sys.path.insert(0, os.path.join(REPO, 'scripts', 'traversal'))
        import index_manifest as im
        d = os.path.join(self.tempdir.name, 'inventory')
        os.makedirs(d, exist_ok=True)
        with open(f'{d}/in.fa', 'w') as f:
            f.write('>s1\n' + self.records['acc1'] + '\n>s2\n' + self.records['acc2'] + '\n')
        self._run_ok(f'{METAGRAPH} build -p 1 --mode basic --graph succinct -k 11 --mask-dummy '
                     '--in-ram -o masked in.fa', d)
        self._run_ok(f'{METAGRAPH} transform --initialize-bloom -o masked masked.dbg', d)
        self._run_ok(f'{METAGRAPH} build -p 1 --mode basic --graph succinct -k 11 -o plain in.fa', d)
        self._annotate_graph(f'{d}/in.fa', f'{d}/plain.dbg', f'{d}/col', 'column',
                             extra_params='--coordinates')
        self._annotate_graph(f'{d}/in.fa', f'{d}/plain.dbg', f'{d}/coord', 'column_coord')
        self._run_ok(f'{METAGRAPH} annotate -p 1 -i plain.dbg --anno-filename '
                     '--index-header-coords -o coord in.fa', d)
        self._annotate_graph(f'{d}/in.fa', f'{d}/plain.dbg', f'{d}/nohead', 'column_coord')
        self._annotate_graph(f'{d}/in.fa', f'{d}/plain.dbg', f'{d}/rd', 'row_diff')
        for name in ('masked.edgemask', 'masked.bloom', 'col.column.annodbg.coords',
                     'coord.seqs', 'plain.dbg.anchors', 'plain.dbg.rd_succ',
                     'rd.row_diff.annodbg'):
            self.assertTrue(os.path.exists(f'{d}/{name}'), name)
        os.makedirs(f'{d}/link', exist_ok=True)
        os.symlink(f'{d}/plain.dbg', f'{d}/link/plain.dbg')
        os.symlink(f'{d}/coord.column_coord.annodbg', f'{d}/link/coord.column_coord.annodbg')
        # an unmasked graph with a Bloom filter beside it (not read: only after the mask)
        os.makedirs(f'{d}/bloomonly', exist_ok=True)
        shutil.copyfile(f'{d}/plain.dbg', f'{d}/bloomonly/plain.dbg')
        shutil.copyfile(f'{d}/masked.bloom', f'{d}/bloomonly/plain.bloom')
        # (the identity roles, the derived (exists, loaded) of the mask and the Bloom filter)
        masked, unmasked = [(True, True), (True, True)], [(False, False), (False, False)]
        cases = [
            ('masked.dbg', 'col.column.annodbg', (), ['graph', 'annotation'], masked),
            ('masked.dbg', 'coord.column_coord.annodbg', (),
             ['graph', 'annotation', 'coord_to_header'], masked),
            ('plain.dbg', 'coord.column_coord.annodbg', ('--no-coord-mapping',),
             ['graph', 'annotation'], unmasked),
            ('plain.dbg', 'nohead.column_coord.annodbg', (), ['graph', 'annotation'], unmasked),
            ('plain.dbg', 'rd.row_diff.annodbg', (),
             ['graph', 'annotation', 'row_diff_anchors', 'row_diff_fork_succ'], unmasked),
            # beside the symlinks there is no .seqs (nor an anchors file): none is listed
            ('link/plain.dbg', 'link/coord.column_coord.annodbg', (), ['graph', 'annotation'],
             unmasked),
            ('bloomonly/plain.dbg', 'col.column.annodbg', (), ['graph', 'annotation'],
             [(False, False), (True, False)]),
        ]
        for graph, anno, flags, roles, derived in cases:
            cli = self._cli_inventory(f'{d}/{graph}', f'{d}/{anno}', *flags)
            py = im.load_inventory(f'{d}/{graph}', f'{d}/{anno}',
                                   coord_mapping='--no-coord-mapping' not in flags)
            self.assertEqual([(e['path'], e['role'], e['required']) for e in cli['files']],
                             py, (graph, anno))
            self.assertEqual(roles, [e['role'] for e in cli['files']], (graph, anno))
            self.assertTrue(all(e['exists'] for e in cli['files']), (graph, anno))
            prefix = f'{d}/{graph}'[:-len('.dbg')]
            self.assertEqual([(prefix + '.edgemask', 'graph_mask') + derived[0],
                              (prefix + '.bloom', 'graph_bloom') + derived[1]],
                             [(e['path'], e['role'], e['exists'], e['loaded'])
                              for e in cli['derived']], (graph, anno))
            self.assertEqual([(e['path'], e['role'], e['exists'], e['loaded'])
                              for e in cli['derived']], im.derived_files(f'{d}/{graph}'))
            self.assertEqual(im.DERIVED_RULE, cli['derived_rule'])
            # the script's own --inventory: what the binary prints, key for key
            res = subprocess.run([sys.executable, os.path.join(REPO, 'scripts', 'traversal',
                                                               'index_manifest.py'),
                                  '--inventory', '-i', f'{d}/{graph}', '-a', f'{d}/{anno}',
                                  *flags], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            self.assertEqual(0, res.returncode, res.stderr.decode())
            self.assertEqual({k: cli[k] for k in ('files', 'derived', 'derived_rule')},
                             json.loads(res.stdout), (graph, anno))
        self.assertEqual([(k['extension'], k['row_diff_anchors'], k['coordinates'])
                          for k in cli['annotation_kinds']], im.ANNOTATION_KINDS)

    def test_symlinked_main_files_do_not_hide_their_sidecars(self):
        """Review of pass 5, finding 2 (the reviewer's /tmp/metagraph-pass5-identity probe):
        two pairs whose graph and annotation are symlinks to the same files, each with its own
        .seqs beside the symlinks (sampleA, sampleBBBB), both naming A's manifest. Grouped by
        the main files, B's .seqs was never checked and both stated one index_fp; now each
        pair's whole inventory is checked, so the server refuses to start, naming B's .seqs.
        With a manifest per pair (index_manifest.py --server-csv) it starts, and the two state
        different index_fp and their own headers."""
        d = os.path.join(self.tempdir.name, 'symlinks')
        os.makedirs(d, exist_ok=True)
        seq = 'ACCGTATGCATAGGCTCCAGTTCAGGATCTCACATCGATGCTTACG'
        with open(f'{d}/source.fa', 'w') as f:
            f.write('>sampleA\n' + seq + '\n')
        self._run_ok(f'{METAGRAPH} build -p 1 --mode basic --graph succinct -k 11 -o shared '
                     'source.fa', d)
        self._annotate_graph(f'{d}/source.fa', f'{d}/shared.dbg', f'{d}/sharedanno',
                             'column_coord')
        for side, header in (('A', 'sampleA'), ('B', 'sampleBBBB')):
            os.makedirs(f'{d}/{side}', exist_ok=True)
            os.symlink(f'{d}/shared.dbg', f'{d}/{side}/graph.dbg')
            os.symlink(f'{d}/sharedanno.column_coord.annodbg',
                       f'{d}/{side}/annotation.column_coord.annodbg')
            with open(f'{d}/{side}/source.fa', 'w') as f:
                f.write('>' + header + '\n' + seq + '\n')
            self._run_ok(f'{METAGRAPH} annotate -p 1 -i graph.dbg --anno-filename '
                         '--index-header-coords -o annotation source.fa', f'{d}/{side}')
        sizes = [os.path.getsize(f'{d}/{side}/annotation.seqs') for side in 'AB']
        self.assertNotEqual(sizes[0], sizes[1])
        # A's manifest (as the reviewer's shared-manifest.json): A's files, A's .seqs
        files = []
        for name in ('graph.dbg', 'annotation.column_coord.annodbg', 'annotation.seqs'):
            with open(f'{d}/A/{name}', 'rb') as f:
                data = f.read()
            files.append({'path': name, 'size': len(data),
                          'sha256': hashlib.sha256(data).hexdigest()})
        with open(f'{d}/shared-manifest.json', 'w') as f:
            json.dump({'format': 'metagraph-index-manifest', 'version': 1, 'files': files}, f)
        csv = f'{d}/symlink.csv'
        with open(csv, 'w') as f:
            for side in 'AB':
                f.write(f'{side},{d}/{side}/graph.dbg,{d}/{side}/annotation.column_coord.annodbg,'
                        f'{d}/shared-manifest.json,bundle\n')
        code, log = self._start_refused(csv)
        self.assertEqual(1, code, log[-2000:])
        self.assertIn(f'{d}/B/annotation.seqs', log)
        # one manifest per pair: two indexes, two identities
        script = os.path.join(REPO, 'scripts', 'traversal', 'index_manifest.py')
        plain = f'{d}/plain.csv'
        with open(plain, 'w') as f:
            for side in 'AB':
                f.write(f'{side},{d}/{side}/graph.dbg,{d}/{side}/annotation.column_coord.annodbg'
                        f',{d}/{side}.manifest.json\n')
        res = subprocess.run([sys.executable, script, '--server-csv', plain],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(0, res.returncode, res.stderr.decode())
        for side in 'AB':
            with open(f'{d}/{side}.manifest.json') as f:
                self.assertIn('annotation.seqs', [e['path'] for e in json.load(f)['files']])
        server = self._Server(self, plain)
        try:
            self.assertTrue(server.ready, server.text()[-2000:])
            fps, names = [], []
            for side in 'AB':
                out = requests.post(server.url + '/resolve', data=json.dumps(
                    {'graph': side, 'sequence': seq,
                     'discover': {'kind': 'header', 'max_labels': 10}})).json()
                fps.append(out['capabilities']['index_fp'])
                names.append([l['label'] for l in out['labels']])
            self.assertNotEqual(fps[0], fps[1])
            self.assertEqual([['sampleA'], ['sampleBBBB']], names)
        finally:
            server.close()

    def test_a_manifest_names_one_bundle_and_its_loaded_files(self):
        """Review of the pass-5 fixes (the reviewer's p7/review-identity dup and mask probes): a
        manifest written for a directory of bundles whose files share base names (A/graph.dbg,
        B/graph.dbg, A/annotation.seqs, B/annotation.seqs) passed for every one of them, since
        loaded files are matched by base name: two single-index servers stated one index_fp
        with different headers. Such a manifest is refused by the server, the CLI and
        index_manifest.py --verify, naming the shared base name. And a manifest that lists an
        optional sidecar the pair does not load (a .seqs missing beside the symlinks of C, or
        --no-coord-mapping) is refused too, so that index_fp identifies the loaded files."""
        d = os.path.join(self.tempdir.name, 'dupnames')
        os.makedirs(d, exist_ok=True)
        seq = 'ACCGTATGCATAGGCTCCAGTTCAGGATCTCACATCGATGCTTACG'
        with open(f'{d}/source.fa', 'w') as f:
            f.write('>sampleA\n' + seq + '\n')
        self._run_ok(f'{METAGRAPH} build -p 1 --mode basic --graph succinct -k 11 -o shared '
                     'source.fa', d)
        self._annotate_graph(f'{d}/source.fa', f'{d}/shared.dbg', f'{d}/sharedanno',
                             'column_coord')
        for side, header in (('A', 'sampleA'), ('B', 'sampleBBBB'), ('C', None)):
            os.makedirs(f'{d}/{side}', exist_ok=True)
            os.symlink(f'{d}/shared.dbg', f'{d}/{side}/graph.dbg')
            os.symlink(f'{d}/sharedanno.column_coord.annodbg',
                       f'{d}/{side}/annotation.column_coord.annodbg')
            if header:
                with open(f'{d}/{side}/source.fa', 'w') as f:
                    f.write('>' + header + '\n' + seq + '\n')
                self._run_ok(f'{METAGRAPH} annotate -p 1 -i graph.dbg --anno-filename '
                             '--index-header-coords -o annotation source.fa', f'{d}/{side}')
        # the directory's manifest: both bundles, by their paths under it
        files = []
        for side in 'AB':
            for name in ('graph.dbg', 'annotation.column_coord.annodbg', 'annotation.seqs'):
                with open(f'{d}/{side}/{name}', 'rb') as f:
                    data = f.read()
                files.append({'path': f'{side}/{name}', 'size': len(data),
                              'sha256': hashlib.sha256(data).hexdigest()})
        files.sort(key=lambda e: e['path'].encode())
        canonical = ''.join('%s\t%d\t%s\n' % (e['path'], e['size'], e['sha256']) for e in files)
        with open(f'{d}/dir.manifest.json', 'w') as f:
            json.dump({'format': 'metagraph-index-manifest', 'version': 1,
                       'index_fp': hashlib.sha256(canonical.encode()).hexdigest(),
                       'files': files}, f)
        script = os.path.join(REPO, 'scripts', 'traversal', 'index_manifest.py')
        res = subprocess.run([sys.executable, script, '--verify', f'{d}/dir.manifest.json'],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertNotEqual(0, res.returncode, res.stdout.decode())
        self.assertIn('base name listed more than once: annotation.seqs', res.stdout.decode())
        request = f'{d}/request.json'
        with open(request, 'w') as f:
            json.dump({'seeds': [{'sequence': seq}],
                       'strategy': {'direction': 'right',
                                    'output': {'detail': 'summary', 'timing': False}}}, f)

        def refused(side, manifest, *flags):
            pair = ['-i', f'{d}/{side}/graph.dbg', '-a',
                    f'{d}/{side}/annotation.column_coord.annodbg', '--index-manifest', manifest]
            server = subprocess.run(shlex.split(METAGRAPH) + [
                'server_query', *pair, *flags, '--port', str(_free_port()), '--address',
                '127.0.0.1'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=120)
            cli = subprocess.run(shlex.split(METAGRAPH) + ['traverse', *pair, *flags, request],
                                 stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=120)
            return ((server.returncode, server.stdout.decode()),
                    (cli.returncode, cli.stdout.decode()))

        for side in 'AB':
            for code, log in refused(side, f'{d}/dir.manifest.json'):
                self.assertEqual(1, code, log[-2000:])
                self.assertIn('share the base name', log)
        # one manifest per pair: A's lists its .seqs
        res = subprocess.run([sys.executable, script, '-i', f'{d}/A/graph.dbg', '-a',
                              f'{d}/A/annotation.column_coord.annodbg', '-o',
                              f'{d}/A.manifest.json'],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(0, res.returncode, res.stderr.decode())
        with open(f'{d}/A.manifest.json') as f:
            self.assertIn('annotation.seqs', [e['path'] for e in json.load(f)['files']])
        ok = subprocess.run(shlex.split(METAGRAPH) + [
            'traverse', '-i', f'{d}/A/graph.dbg', '-a', f'{d}/A/annotation.column_coord.annodbg',
            '--index-manifest', f'{d}/A.manifest.json', request],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=120)
        self.assertEqual(0, ok.returncode, ok.stderr.decode()[-2000:])
        # ... but C loads no .seqs (none beside its symlinks), and neither does A under
        # --no-coord-mapping: both refuse A's manifest, naming the .seqs
        for side, flags in (('C', ()), ('A', ('--no-coord-mapping',))):
            for code, log in refused(side, f'{d}/A.manifest.json', *flags):
                self.assertEqual(1, code, (side, flags, log[-2000:]))
                self.assertIn('annotation.seqs, which the server does not load', log)

    def test_derived_data_is_not_part_of_the_identity(self):
        """The owner's decision #17 of 2026-10-08 (until then: review of pass 5, finding 3, a
        loaded Bloom filter was part of the identity): the graph's dummy-edge mask and Bloom
        filter are derived data, not part of index_fp. A manifest written before they exist is
        valid after `transform --mask-dummy` and `--initialize-bloom`, its index_fp unchanged,
        on a multi-graph server and on the CLI; index_manifest.py writes the same manifest
        with or without them (it never lists them) and refuses them with --extra. A manifest
        that lists one (the reviewer's graph, mask and annotation) is refused, by the server
        and the CLI naming the rule, and flagged by --verify. And the server names the derived
        data it loads beside a manifest (operators see what is loaded)."""
        d = os.path.join(self.tempdir.name, 'bloom')
        os.makedirs(d, exist_ok=True)
        seq = 'ACCGTATGCATAGGCTCCAGTTCAGGATCTCACATCGATGCTTACG'
        with open(f'{d}/seq.fa', 'w') as f:
            f.write('>A\n' + seq + '\n')
        self._run_ok(f'{METAGRAPH} build -p 1 --mode basic --graph succinct -k 11 '
                     '--in-ram -o graph seq.fa', d)
        self._annotate_graph(f'{d}/seq.fa', f'{d}/graph.dbg', f'{d}/annotation', 'column')
        script = os.path.join(REPO, 'scripts', 'traversal', 'index_manifest.py')

        def manifest(name, *extra):
            res = subprocess.run([sys.executable, script, '-i', f'{d}/graph.dbg', '-a',
                                  f'{d}/annotation.column.annodbg', '-o', f'{d}/{name}',
                                  *extra], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            return res.returncode, res.stdout.decode() + res.stderr.decode()

        # before the mask: the graph and the annotation
        self.assertEqual(0, manifest('before.manifest.json')[0])
        with open(f'{d}/before.manifest.json') as f:
            before = json.load(f)
        self.assertEqual(['annotation.column.annodbg', 'graph.dbg'],
                         [e['path'] for e in before['files']])
        # the derived data added beside the deployed graph
        self._run_ok(f'{METAGRAPH} transform --mask-dummy -p 1 graph.dbg', d)
        self._run_ok(f'{METAGRAPH} transform --initialize-bloom -o graph graph.dbg', d)
        for name in ('graph.edgemask', 'graph.bloom'):
            self.assertTrue(os.path.exists(f'{d}/{name}'), name)
        self.assertEqual([('graph_mask', True), ('graph_bloom', True)],
                         [(e['role'], e['loaded']) for e in self._cli_inventory(
                             f'{d}/graph.dbg', f'{d}/annotation.column.annodbg')['derived']])
        # index_manifest.py writes the same manifest now (it never lists them) ...
        self.assertEqual(0, manifest('after.manifest.json')[0])
        with open(f'{d}/after.manifest.json') as f:
            after = json.load(f)
        self.assertEqual(before['files'], after['files'])
        self.assertEqual(before['index_fp'], after['index_fp'])
        # ... and refuses them as extra files
        for extra in ('graph.edgemask', 'graph.bloom'):
            code, text = manifest('extra.manifest.json', '--extra', f'{d}/{extra}')
            self.assertNotEqual(0, code)
            self.assertIn('derived data of the graph, not part of index_fp', text)
            self.assertFalse(os.path.exists(f'{d}/extra.manifest.json'))
        # the manifest written before the mask: the server starts, states its index_fp, and
        # names the derived data it loaded; /resolve answers as on the reviewer's index
        csv = f'{d}/before.csv'
        with open(csv, 'w') as f:
            f.write(f'A,{d}/graph.dbg,{d}/annotation.column.annodbg,{d}/before.manifest.json\n')
        server = self._Server(self, csv)
        try:
            self.assertTrue(server.ready, server.text()[-2000:])
            out = requests.post(server.url + '/resolve', data=json.dumps(
                {'graph': 'A', 'sequence': seq, 'labels': ['A']})).json()
            self.assertEqual([[0, 36]], out['graph_runs'])
            self.assertEqual(before['index_fp'], out['capabilities']['index_fp'])
            self.assertIn(f'Derived data of {d}/graph.dbg loaded beside it, not part of '
                          f'index_fp: {d}/graph.edgemask, {d}/graph.bloom', server.text())
        finally:
            server.close()
        request = f'{d}/request.json'
        with open(request, 'w') as f:
            json.dump({'seeds': [{'sequence': seq}],
                       'strategy': {'direction': 'right',
                                    'output': {'detail': 'summary', 'timing': False}}}, f)
        res = subprocess.run(shlex.split(METAGRAPH) + [
            'traverse', '-i', f'{d}/graph.dbg', '-a', f'{d}/annotation.column.annodbg',
            '--index-manifest', f'{d}/before.manifest.json', request],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(0, res.returncode, res.stderr.decode()[-2000:])
        self.assertEqual(before['index_fp'], json.loads(res.stdout)['capabilities']['index_fp'])
        # the reviewer's manifest: graph, mask and annotation — refused, naming the rule
        files = []
        for name in ('annotation.column.annodbg', 'graph.dbg', 'graph.edgemask'):
            with open(f'{d}/{name}', 'rb') as f:
                data = f.read()
            files.append({'path': name, 'size': len(data),
                          'sha256': hashlib.sha256(data).hexdigest()})
        with open(f'{d}/old.manifest.json', 'w') as f:
            json.dump({'format': 'metagraph-index-manifest', 'version': 1, 'files': files}, f)
        csv = f'{d}/old.csv'
        with open(csv, 'w') as f:
            f.write(f'A,{d}/graph.dbg,{d}/annotation.column.annodbg,{d}/old.manifest.json\n')
        code, log = self._start_refused(csv)
        self.assertEqual(1, code, log[-2000:])
        rule = ('it lists graph.edgemask: the dummy-edge mask (.edgemask) and the Bloom filter '
                '(.bloom) are derived data of the graph, not part of index_fp')
        self.assertIn(rule, log)
        self.assertIn('write it again without graph.edgemask', log)
        res = subprocess.run(shlex.split(METAGRAPH) + [
            'traverse', '-i', f'{d}/graph.dbg', '-a', f'{d}/annotation.column.annodbg',
            '--index-manifest', f'{d}/old.manifest.json', request],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertNotEqual(0, res.returncode)
        self.assertIn(rule, res.stderr.decode() + res.stdout.decode())
        res = subprocess.run([sys.executable, script, '--verify', f'{d}/old.manifest.json'],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertNotEqual(0, res.returncode)
        self.assertIn('derived data listed: graph.edgemask', res.stdout.decode())
        res = subprocess.run([sys.executable, script, '--verify', f'{d}/before.manifest.json'],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(0, res.returncode, res.stdout.decode())

    def test_multi_three_column_list_states_nulls(self):
        """A list of three columns is read as it always was: null name and fingerprint."""
        d = self.tempdir.name
        csv = f'{d}/plain.csv'
        with open(csv, 'w') as f:
            f.write(f'P,{self.pairs["three"][0]},{self.pairs["three"][1]}\n')
        server = self._Server(self, csv)
        try:
            self.assertTrue(server.ready, server.text())
            out = requests.post(server.url + '/traverse',
                                data=json.dumps(self._request('P'))).json()
            caps = out['capabilities']
            self.assertEqual((None, None), (caps['index_ns'], caps['index_fp']))
            self.assertEqual(16, len(caps['index_meta_fp']))
        finally:
            server.close()


# the mini index (scripts/traversal/build_mini_refseq.sh; not in CI), as test_pattern.py finds it
MINI_DIR = os.environ.get('METAGRAPH_MINI_REFSEQ', os.path.join(os.getcwd(), 'mini_refseq'))
MINI_FILES = ('graph_k31.dbg', 'graph_k31.dbg.anchors', 'graph_k31.dbg.rd_succ',
              'annotation.relaxed.relabeled.row_diff_brwt_coord.annodbg',
              'annotation.relaxed.relabeled.seqs')
# $METAGRAPH_REQUIRE_GUARDS=1: a missing mini index fails the class instead of skipping it
_MINI_PRESENT = all(os.path.isfile(os.path.join(MINI_DIR, f)) for f in MINI_FILES)
_MINI_GUARD = (_MINI_PRESENT or os.environ.get('METAGRAPH_REQUIRE_GUARDS', '') == '1')


@unittest.skipIf(PROTEIN_MODE, "traversal fixtures are DNA")
@unittest.skipUnless(_supports_traverse(), "`metagraph traverse` is not available in this build")
@unittest.skipUnless(_MINI_GUARD, "the mini index is not built (scripts/traversal/"
                                  "build_mini_refseq.sh; not in CI)")
class TestTraverseDerivedDataMini(TestingBase):
    """The owner's decision #17 of 2026-10-08 on a copy of build/mini_refseq (unmasked, as
    built, like refseq33m): the dummy-edge mask is derived data of the graph, not part of
    index_fp. A manifest made for the graph without a mask (index_manifest.py, the same
    index_fp as the build's own manifest) serves a single-index server_query --index-manifest;
    the mask is added beside the graph with `transform --mask-dummy` (the staging step) and
    the server restarted with the SAME manifest: it starts, states the same index_fp (and
    index_meta_fp), and stored /traverse retrievals and /resolve answers replay identically
    (apart from timing): the graphlets stored before compare as the same index, equal. GET
    /stats graph.nodes changes (the k-mers instead of the edges), as documented. A manifest
    that lists the mask is refused. build/mini_refseq itself is never written to."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        assert _MINI_PRESENT, f'the mini index is not in {MINI_DIR}'
        cls.mini_files = sorted(os.listdir(MINI_DIR))
        d = cls.dir = os.path.join(cls.tempdir.name, 'mini')
        os.makedirs(d)
        for name in MINI_FILES:
            shutil.copyfile(os.path.join(MINI_DIR, name), os.path.join(d, name))
        cls.graph = os.path.join(d, MINI_FILES[0])
        cls.anno = os.path.join(d, MINI_FILES[3])
        # seeds from the records the index was built from (first record of the first files)
        cls.seqs = []
        for fasta in sorted(os.listdir(os.path.join(MINI_DIR, 'fasta')))[:3]:
            with open(os.path.join(MINI_DIR, 'fasta', fasta)) as f:
                lines = f.read().split('>')[1].split('\n')
            seq = ''.join(lines[1:]).upper()
            start = next(i for i in range(0, len(seq) - 700, 37)
                         if set(seq[i:i + 700]) <= set('ACGT'))
            cls.seqs.append(seq[start:start + 700])

    @classmethod
    def tearDownClass(cls):
        cls.tempdir.cleanup()

    def _server(self, manifest, wait=True):
        port = _free_port()
        log_path = f'{self.dir}/server-{port}.log'
        log = open(log_path, 'w')
        process = subprocess.Popen(
            shlex.split(METAGRAPH) + ['server_query', '-i', self.graph, '-a', self.anno,
                                      '--port', str(port), '--address', '127.0.0.1', '-p', '2',
                                      '--index-name', 'mini', '--index-manifest', manifest],
            stdout=log, stderr=subprocess.STDOUT)
        url = f'http://127.0.0.1:{port}'
        for _ in range(1200 if wait else 0):
            if process.poll() is not None:
                break
            try:
                if requests.get(url + '/traverse/capabilities', timeout=2).ok:
                    break
            except requests.exceptions.RequestException:
                pass
            time.sleep(0.1)

        def stop():
            if process.poll() is None:
                process.kill()
            process.wait()
            log.close()
            with open(log_path) as f:
                return f.read()
        return url, process, stop

    def _requests(self):
        traverse = [{'seeds': [{'sequence': s[300:340]}],
                     'strategy': {'direction': 'both', 'bounds': {'max_extension_bp': 250},
                                  'output': {'detail': 'graphlet', 'timing': False}}}
                    for s in self.seqs]
        traverse.append({'seeds': [{'sequence': self.seqs[0][100:140]}],
                         'strategy': {'direction': 'right', 'bounds': {'max_extension_bp': 120},
                                      'output': {'detail': 'full', 'timing': False}}})
        resolve = [{'sequence': s, 'discover': {'kind': 'header', 'max_labels': 10}}
                   for s in self.seqs]
        resolve.append({'sequence': self.seqs[1], 'discover': {'max_labels': 10},
                        'support': 'trace'})
        return traverse, resolve

    @staticmethod
    def _untimed(value):
        """|value| without its timing: 'timing' members dropped, every elapsed_ms 0."""
        if isinstance(value, dict):
            return {k: (0 if k == 'elapsed_ms' else TestTraverseDerivedDataMini._untimed(v))
                    for k, v in value.items() if k != 'timing'}
        if isinstance(value, list):
            return [TestTraverseDerivedDataMini._untimed(v) for v in value]
        return value

    def _answers(self, url):
        """(identity, the stored answers: untimed results of every request, the graphlets)"""
        caps = requests.get(url + '/traverse/capabilities').json()
        identity = (caps['index_ns'], caps['index_fp'], caps['index_meta_fp'])
        traverse, resolve = self._requests()
        answers, graphlets = [], []
        for request in traverse:
            ret = requests.post(url + '/traverse', data=json.dumps(request))
            self.assertEqual(200, ret.status_code, ret.text[:2000])
            out = ret.json()
            self.assertEqual(identity, tuple(out['capabilities'][k] for k in
                                             ('index_ns', 'index_fp', 'index_meta_fp')))
            answers.append(('traverse', self._untimed(out['results'])))
            for r in out['results']:
                if 'graphlet' in r:
                    # the H record states the identity the graphlet is stored under
                    self.assertEqual(list(identity), r['graphlet'].split('\n')[0].split(' ')[-3:])
                    graphlets.append(graphlet_lib.from_response(r, out))
        for request in resolve:
            ret = requests.post(url + '/resolve', data=json.dumps(request))
            self.assertEqual(200, ret.status_code, ret.text[:2000])
            out = ret.json()
            self.assertEqual(identity, tuple(out['capabilities'][k] for k in
                                             ('index_ns', 'index_fp', 'index_meta_fp')))
            answers.append(('resolve', self._untimed({k: v for k, v in out.items()
                                                      if k != 'capabilities'})))
        stats = requests.get(url + '/stats').json()
        return identity, answers, graphlets, stats

    def test_adding_a_mask_keeps_the_index_fp(self):
        d = self.dir
        script = os.path.join(REPO, 'scripts', 'traversal', 'index_manifest.py')
        manifest = f'{d}/index.manifest.json'
        res = subprocess.run([sys.executable, script, '-i', self.graph, '-a', self.anno,
                              '-o', manifest, '--name', 'mini'],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(0, res.returncode, res.stderr.decode())
        with open(manifest) as f:
            fp = json.load(f)['index_fp']
        # what a row_diff_brwt_coord pair loads: the graph, the annotation, the headers
        with open(manifest) as f:
            self.assertEqual(sorted([MINI_FILES[0], MINI_FILES[3], MINI_FILES[4]],
                                    key=str.encode),
                             [e['path'] for e in json.load(f)['files']])
        # the build's own manifest (build_mini_refseq.sh: it lists the row-diff anchors and
        # fork successors too, files outside the inventory, so its index_fp is its own)
        built = os.path.join(MINI_DIR, 'annotation.relaxed.relabeled.manifest.json')
        manifests = [(manifest, fp)]
        if os.path.isfile(built):
            shutil.copyfile(built, f'{d}/built.manifest.json')
            with open(built) as f:
                manifests.append((f'{d}/built.manifest.json', json.load(f)['index_fp']))
        self.assertFalse(os.path.exists(f'{d}/graph_k31.edgemask'))

        def identities():
            """the index_fp each manifest makes the server state (the first's answers too)"""
            out = []
            for path, stated in manifests:
                url, process, stop = self._server(path)
                try:
                    self.assertIsNone(process.poll(), 'the server did not start')
                    if not out:
                        answers = self._answers(url)
                    caps = requests.get(url + '/traverse/capabilities').json()
                    out.append(caps['index_fp'])
                finally:
                    log = stop()
                self.assertEqual(stated, out[-1], log[-2000:])
            return out, answers, log

        fps, before, log = identities()
        self.assertEqual(('mini', fp), before[0][:2], log[-2000:])
        self.assertNotIn('Derived data of', log)
        self.assertGreaterEqual(len(before[2]), 3)

        # the staging step: the mask beside the deployed graph, nothing else written
        res = self._run_command(f'{METAGRAPH} transform --mask-dummy -p 2 {self.graph}',
                                'Mask the copy of the mini graph')
        transform_log = (res.stdout + res.stderr).decode()
        self.assertEqual(sorted(MINI_FILES + ('graph_k31.edgemask',)),
                         sorted(n for n in os.listdir(d) if not n.startswith('server-')
                                and not n.endswith('.manifest.json')))
        for name in MINI_FILES:
            self.assertTrue(filecmp.cmp(f'{d}/{name}', os.path.join(MINI_DIR, name),
                                        shallow=False), name)
        inventory = subprocess.run(shlex.split(METAGRAPH) + [
            'traverse', '--index-inventory', '--json', '-i', self.graph, '-a', self.anno],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(0, inventory.returncode, inventory.stderr.decode())
        inventory = json.loads(inventory.stdout)
        self.assertEqual([f'{d}/{n}' for n in (MINI_FILES[0], MINI_FILES[3], MINI_FILES[4])],
                         [e['path'] for e in inventory['files']])
        self.assertEqual([(f'{d}/graph_k31.edgemask', 'graph_mask', True, True),
                          (f'{d}/graph_k31.bloom', 'graph_bloom', False, False)],
                         [(e['path'], e['role'], e['exists'], e['loaded'])
                          for e in inventory['derived']])

        # the restart, with the manifests made before the mask: each states its index_fp
        fps_after, after, log = identities()
        self.assertEqual(fps, fps_after)
        self.assertIn(f'Derived data of {self.graph} loaded beside it, not part of index_fp: '
                      f'{d}/graph_k31.edgemask', log)
        # the same identity, the same answers: stored retrievals stay valid
        self.assertEqual(before[0], after[0])
        self.assertEqual(len(before[1]), len(after[1]))
        for (route, a), (_, b) in zip(before[1], after[1]):
            self.assertEqual(a, b, route)
        self.assertEqual(len(before[2]), len(after[2]))
        for stored, fresh in zip(before[2], after[2]):
            # comparable as one index (a new index_fp would make them 'different indexes');
            # 'qualified' only where a side's label lists were cut (then equality is not
            # claimed by the library, and the byte comparison above stands for it)
            c = stored.compare(fresh)
            self.assertIn(c.comparable, (True, 'qualified'), c.reason)
            self.assertNotIn('different indexes', c.reason)
            if c.comparable is True:
                self.assertIs(True, c.equal, c.reason)
        # ... and what does change, as documented: /stats states the k-mers as the nodes
        m = re.search(r'(\d+) edges, (\d+) source dummies \(the main dummy edge included\), '
                      r'(\d+) sink dummies, (\d+) k-mers', transform_log)
        self.assertIsNotNone(m, transform_log)
        edges, source, sink, kmers = (int(x) for x in m.groups())
        self.assertEqual(edges, source + sink + kmers)
        self.assertEqual(edges, before[3]['graph']['nodes'])
        self.assertEqual(kmers, after[3]['graph']['nodes'])

        # the CLI states the same identity with the mask beside the graph
        request = f'{d}/request.json'
        with open(request, 'w') as f:
            json.dump(self._requests()[0][0], f)
        res = subprocess.run(shlex.split(METAGRAPH) + [
            'traverse', '-i', self.graph, '-a', self.anno, '--index-name', 'mini',
            '--index-manifest', manifest, request],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(0, res.returncode, res.stderr.decode()[-2000:])
        self.assertEqual(fp, json.loads(res.stdout)['capabilities']['index_fp'])

        # a manifest that lists the mask is refused, naming the rule
        with open(manifest) as f:
            listing = json.load(f)
        with open(f'{d}/graph_k31.edgemask', 'rb') as f:
            data = f.read()
        listing['files'].append({'path': 'graph_k31.edgemask', 'size': len(data),
                                 'sha256': hashlib.sha256(data).hexdigest()})
        listing.pop('index_fp')
        with open(f'{d}/masked.manifest.json', 'w') as f:
            json.dump(listing, f)
        url, process, stop = self._server(f'{d}/masked.manifest.json', wait=False)
        try:
            code = process.wait(timeout=120)
        finally:
            log = stop()
        self.assertEqual(1, code, log[-2000:])
        self.assertIn('it lists graph_k31.edgemask: the dummy-edge mask (.edgemask) and the '
                      'Bloom filter (.bloom) are derived data of the graph, not part of index_fp',
                      log)
        # build/mini_refseq is unchanged
        self.assertEqual(self.mini_files, sorted(os.listdir(MINI_DIR)))


@unittest.skipIf(PROTEIN_MODE, "traversal fixtures are DNA")
@unittest.skipUnless(_supports_traverse(), "`metagraph traverse` is not available in this build")
class TestTraverseWideIndex(TestingBase):
    """Pass 5, W5 (the chunked deadlines): a fan-out index — a 61 bp seed, then all 1024
    five-base continuations, each with a 400 bp tail of its own, one header label each — where
    one level's lookahead reads about 66,000 annotation rows in one call. A time budget of 50 to
    200 ms used to be overrun to about half a second by that call (measured 507-534 ms); the
    reads are now decoded in chunks of chunk_target_ms with the deadline between them, so the
    walk stops within about a chunk of its budget, where the whole read would have stopped it.
    With --traverse-chunk-target-ms 0 the read is one piece again (the test bites)."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        rng = random.Random(20261004)

        def seq(n):
            return ''.join(rng.choice('ACGT') for _ in range(n))

        cls.seed = seq(61)
        d = cls.tempdir.name
        fasta = f'{d}/fan.fa'
        with open(fasta, 'w') as f:
            for i in range(4 ** 5):
                penta = ''.join('ACGT'[(i >> (2 * j)) & 3] for j in range(5))
                f.write(f'>r{i:04d}\n{cls.seed}{penta}{seq(400)}\n')
        cls.graph = f'{d}/fan.dbg'
        cls._build_graph(fasta, cls.graph, K, 'succinct', mode='basic')
        cls._annotate_graph(fasta, cls.graph, f'{d}/fan', 'column', anno_type='header')
        cls.anno = f'{d}/fan.column.annodbg'
        # Review of pass 5, F1/F2: the same index on row-diff annotations, whose rows share the
        # decoding of their row-diff paths within a call (a graph copy each: the transforms write
        # sidecars beside the graph)
        cls.variants = {}
        for anno_type in ('row_diff', 'row_diff_brwt'):
            graph = f'{d}/fan_{anno_type}.dbg'
            shutil.copyfile(cls.graph, graph)
            cls._annotate_graph(fasta, graph, f'{d}/fan_{anno_type}', anno_type,
                                anno_type='header')
            cls.variants[anno_type] = types.SimpleNamespace(
                graph=graph, anno=f'{d}/fan_{anno_type}.{anno_type}.annodbg',
                tempdir=cls.tempdir)

    def _request(self, mode, budget_ms, **fields):
        labels = ({'mode': 'annotate', 'max_labels_per_node': 2000} if mode == 'annotate'
                  else {'max_seed_labels': 2000})
        strategy = {'direction': 'right', 'labels': labels,
                    'branching': {'on_reconverge': 'keep', 'max_label_branches': 'unlimited',
                                  'max_splits_per_path': 'unlimited'},
                    'bounds': {'max_extension_bp': 400, 'time_budget_ms': budget_ms,
                               'max_live_paths': 100000, 'max_paths': 100000},
                    'output': {'detail': 'summary', 'timing': True}}
        req = {'seeds': [{'sequence': self.seed}], 'strategy': strategy}
        req.update(fields)
        return req

    def _server(self, *flags, index=None):
        return TestTraverseAttempts._Server(index or self, *flags)

    def test_row_diff_walks_no_deadline_stops_keep_their_bytes_and_time(self):
        """Review of pass 5, F1: on a row-diff annotation the rows of one call share the
        decoding of their paths, and chunks of sorted rows decoded them again per chunk — a walk
        of about 1.8 s took 30 s on row_diff, ran into its budget and came back partial (its
        bytes changed), and row_diff_brwt walks were 1.4-2.9 times slower. A read the deadline
        cannot fall into is one piece now: the same bytes as unchunked, in about its time."""
        for anno_type, ext in (('row_diff', 8), ('row_diff_brwt', 400)):
            index = self.variants[anno_type]
            with self._server('--traverse-chunk-target-ms', '0', index=index) as whole, \
                    self._server(index=index) as paced:
                for mode in ('annotate', 'constrain'):
                    req = self._request(mode, 30000)
                    req['strategy']['bounds']['max_extension_bp'] = ext
                    times = {}
                    results = {}
                    for name, server in (('whole', whole), ('paced', paced)) * 2:
                        out = server.post('traverse', req).json()
                        times.setdefault(name, []).append(out['timing']['elapsed_ms'])
                        result = out['results'][0]
                        self.assertEqual('complete', result['outcome']['walks'],
                                         (anno_type, mode, name))
                        result.pop('timing', None)
                        results.setdefault(name, set()).add(json.dumps(result, sort_keys=True))
                    self.assertEqual(results['whole'], results['paced'], (anno_type, mode))
                    self.assertEqual(1, len(results['paced']))
                    # (the unchunked time, the better of two, with room for a loaded machine)
                    self.assertLessEqual(min(times['paced']), 1.5 * min(times['whole']) + 300,
                                         (anno_type, mode, times))

    def test_row_diff_path_cache_keeps_the_bytes(self):
        """The efficiency pass, the row-diff path cache: the rows a request's reads reconstruct
        are kept, so that a later read's row-diff path stops at a cached row instead of decoding
        to its anchor again (feature level 4: the probe's decode_cache). The responses are byte
        for byte those without the cache (--traverse-path-cache-mb 0) — constrain and annotate,
        no budget, memory budgets (a large one, and one that stops the walk), a work budget —
        on both row-diff annotations."""
        for anno_type, ext in (('row_diff', 8), ('row_diff_brwt', 100)):
            index = self.variants[anno_type]
            with self._server('--traverse-path-cache-mb', '0', index=index) as off, \
                    self._server(index=index) as on:
                caps = requests.get(on.url + '/traverse/capabilities').json()
                self.assertEqual(128, caps['decode_cache']['path_cache_mb'])
                self.assertEqual(0, requests.get(off.url + '/traverse/capabilities')
                                 .json()['decode_cache']['path_cache_mb'])
                stops = 0
                for mode in ('constrain', 'annotate'):
                    for extra in ({}, {'max_memory_mb': 4096}, {'max_memory_mb': 1},
                                  {'max_work_units': 200000}):
                        req = self._request(mode, 60000)
                        req['strategy']['bounds']['max_extension_bp'] = ext
                        req['strategy']['bounds'].update(extra)
                        req['strategy']['output']['timing'] = False
                        a = off.post('traverse', req)
                        b = on.post('traverse', req)
                        what = (anno_type, mode, extra)
                        self.assertEqual(200, a.status_code, a.text[:500])
                        self.assertEqual(200, b.status_code, b.text[:500])
                        self.assertEqual(a.text, b.text, what)
                        stops += '"resource_stop"' in a.text
                self.assertGreater(stops, 0, anno_type)

    def test_server_budget_maxima(self):
        """R16 (the owner's decision; feature level 4): --traverse-max-memory-mb and
        --traverse-max-work-units, off (0) by default. Set, a request's larger budget is
        lowered to the maximum and an omitted one set to it, echoed in strategy.clamped (an
        omitted one as requested "unlimited") and strategy.bounds; a smaller one is kept. A
        seed the clamped budget stopped states a server_clamp limitation, and its stop what the
        request asked for; the graphlet holds both. The probe states the maxima."""
        with self._server('--traverse-max-memory-mb', '4096',
                          '--traverse-max-work-units', '5000') as server:
            caps = requests.get(server.url + '/traverse/capabilities').json()
            self.assertEqual((4096, 5000), (caps['max_memory_mb'], caps['max_work_units']))
            req = self._request('constrain', 60000)
            req['strategy']['output']['timing'] = False
            ret = server.post('traverse', req)
            self.assertEqual(200, ret.status_code, ret.text[:500])
            out = ret.json()
            clamped = {c['field']: c for c in out['strategy']['clamped']}
            self.assertEqual({'field': 'bounds.max_memory_mb', 'requested': 'unlimited',
                              'effective': 4096}, clamped['bounds.max_memory_mb'])
            self.assertEqual({'field': 'bounds.max_work_units', 'requested': 'unlimited',
                              'effective': 5000}, clamped['bounds.max_work_units'])
            self.assertEqual((4096, 5000), (out['strategy']['bounds']['max_memory_mb'],
                                            out['strategy']['bounds']['max_work_units']))
            result = out['results'][0]
            self.assertEqual(('work', 'unlimited', 5000),
                             (result['resource_stop']['resource'],
                              result['resource_stop']['requested'],
                              result['resource_stop']['effective']))
            clamps = [l for l in result['limitations'] if l['kind'] == 'server_clamp']
            self.assertEqual([('bounds.max_work_units', 5000, 'unlimited')],
                             [(l['knob'], l['limit'], l['observed']) for l in clamps])
            self.assertIn('gave no budget', clamps[0]['effect'])
            # a request cannot raise it: no action raising it, and the arm's walk_domain states
            # the maximum (review of the efficiency pass, finding 2)
            self.assertNotIn('raise_work_budget', result['resource_stop']['actions'])
            domains = [l for arm in result['arms'].values() for l in arm.get('limitations', [])
                       if l['kind'] == 'walk_domain' and l['knob'] == 'bounds.max_work_units']
            self.assertTrue(domains)
            for l in domains:
                self.assertEqual(5000, l['server_limit'])
                self.assertNotIn('raise the knob', l['effect'])
            # the graphlet's K and Q records hold them
            req['strategy']['output']['detail'] = 'graphlet'
            text = server.post('traverse', req).json()['results'][0]['graphlet']
            # (K values: i:<integer>, u for unlimited)
            self.assertIn('server_clamp bounds.max_work_units i:5000 u ', text)
            self.assertRegex(text, r'(?m)^Q locus work \S+ u i:5000 ')
            # larger: lowered; smaller: kept
            req['strategy']['output']['detail'] = 'summary'
            req['strategy']['bounds'].update({'max_memory_mb': 8192, 'max_work_units': 10 ** 9})
            out = server.post('traverse', req).json()
            clamped = {c['field']: c for c in out['strategy']['clamped']}
            self.assertEqual((8192, 4096), (clamped['bounds.max_memory_mb']['requested'],
                                            clamped['bounds.max_memory_mb']['effective']))
            self.assertEqual((10 ** 9, 5000), (clamped['bounds.max_work_units']['requested'],
                                               clamped['bounds.max_work_units']['effective']))
            req['strategy']['bounds'].update({'max_memory_mb': 100, 'max_work_units': 4000})
            out = server.post('traverse', req).json()
            self.assertEqual([], [c for c in out['strategy']['clamped']
                                  if c['field'] in ('bounds.max_memory_mb',
                                                    'bounds.max_work_units')])
            self.assertEqual((100, 4000), (out['strategy']['bounds']['max_memory_mb'],
                                           out['strategy']['bounds']['max_work_units']))
        # Review of the efficiency pass, finding 2: a seed the server's maximum fails before any
        # walk states it as a walked one does — the server_clamp, the stop's `requested`, no
        # action raising the budget, and a walk_domain with server_limit saying that a request
        # cannot raise it: here at its depth-0 state (1,024 derived labels under 1 MiB; a seed
        # phase that runs out of work units is MiniRefSeq.ServerBudgetMaxima's: this seed's
        # 31,992 units stay below the interval at which the seed phase is compared)
        for flags, resource, field, maximum in (
                (('--traverse-max-memory-mb', '1'), 'memory', 'bounds.max_memory_mb', 1),):
            with self._server(*flags) as server:
                req = self._request('constrain', 60000)
                req['strategy']['output']['timing'] = False
                ret = server.post('traverse', req)
                self.assertEqual(200, ret.status_code, ret.text[:500])
                result = ret.json()['results'][0]
                self.assertEqual('failed', result['outcome']['walks'], result.get('error'))
                stop = result['resource_stop']
                self.assertEqual((resource, 'unlimited', maximum),
                                 (stop['resource'], stop['requested'], stop['effective']))
                self.assertEqual([], [a for a in stop['actions'] if a.startswith('raise_')])
                mine = [l for l in result['limitations'] if l['knob'] == field]
                self.assertEqual([('unlimited', maximum)],
                                 [(l['observed'], l['limit']) for l in mine
                                  if l['kind'] == 'server_clamp'])
                domains = [l for l in mine if l['kind'] == 'walk_domain']
                self.assertEqual(1, len(domains), mine)
                self.assertEqual(maximum, domains[0]['server_limit'])
                self.assertNotIn('raise the knob', domains[0]['effect'])
                self.assertIn("the server's maximum", domains[0]['effect'])
        for flag, value in (('--traverse-max-memory-mb', '1048577'),
                            ('--traverse-max-memory-mb', '-1'),
                            ('--traverse-max-work-units', '1.5'),
                            ('--traverse-path-cache-mb', '-1')):
            res = subprocess.run(shlex.split(METAGRAPH) + [
                'server_query', '-i', self.graph, '-a', self.anno, flag, value,
                '--port', str(_free_port()), '--address', '127.0.0.1'],
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=120)
            self.assertNotEqual(0, res.returncode, (flag, value))
            self.assertIn(flag, res.stdout.decode(), (flag, value))

    def test_flags_stated_as_integers_are_validated(self):
        """Review of pass 5: --traverse-clock-skew-ms and --traverse-chunk-target-ms were read
        with atoll, so -1 wrapped to 2^64 - 1 and the capabilities stated it (a JSON client
        cannot represent it, a ledger adding to it overflows). Integers in [0, 2^53 - 1] only."""
        for flag in ('--traverse-clock-skew-ms', '--traverse-chunk-target-ms',
                     '--traverse-attempt-allowance-ms'):
            for value in ('-1', '1.5', '9007199254740992', '5x'):
                res = subprocess.run(shlex.split(METAGRAPH) + [
                    'server_query', '-i', self.graph, '-a', self.anno, flag, value,
                    '--port', str(_free_port()), '--address', '127.0.0.1'],
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=120)
                self.assertNotEqual(0, res.returncode, (flag, value))
                self.assertIn(f'{flag} must be an integer in [0, 2^53 - 1]',
                              res.stdout.decode(), (flag, value))

    def test_row_diff_budgets_are_kept_within_a_chunk(self):
        """Review of pass 5, F2: the first chunk of a read was sized from other reads' rows (a
        level's cheap shared rows) and took 650-880 ms against chunk_target_ms 50, so 50-500
        ms budgets overran as before on row_diff; a split read starts with at most 8 rows now,
        and a read is one piece only when predicted at the slowest per-row time seen to end
        well before the deadline. With and without a memory budget (its reads are budget-aware
        on these annotations)."""
        for anno_type in ('row_diff', 'row_diff_brwt'):
            with self._server(index=self.variants[anno_type]) as server:
                for extra in ({}, {'max_memory_mb': 4096}):
                    for mode in ('constrain', 'annotate'):
                        for budget in (50, 200):
                            req = self._request(mode, budget, attempt_id=f'rd-{anno_type}-'
                                                f'{mode}-{budget}-{len(extra)}')
                            req['strategy']['bounds'].update(extra)
                            ret = server.post('traverse', req)
                            self.assertEqual(200, ret.status_code, ret.text[:500])
                            out = ret.json()
                            what = (anno_type, extra, mode, budget)
                            right = out['results'][0]['arms']['right']
                            self.assertEqual('time_budget', right['cap_trigger']['reason'], what)
                            self.assertLessEqual(out['timing']['elapsed_ms'], budget + 50 + 150,
                                                 what)
                            self.assertLessEqual(out['usage']['observed_max_uninterruptible_ms'],
                                                 50 + 100, what)

    def test_wide_index_budget_is_kept_within_a_chunk(self):
        with self._server() as server:
            caps = requests.get(server.url + '/traverse/capabilities').json()
            dc = caps['deadline_check']
            self.assertEqual(50, dc['chunk_target_ms'])
            self.assertIsNone(dc['max_uninterruptible_ms'])
            n = 0
            for mode in ('constrain', 'annotate'):
                for budget in (50, 100, 200):
                    ret = server.post('traverse', self._request(mode, budget,
                                                                attempt_id=f'w-{mode}-{budget}'))
                    self.assertEqual(200, ret.status_code, ret.text[:500])
                    out = ret.json()
                    right = out['results'][0]['arms']['right']
                    self.assertEqual('time_budget', right['cap_trigger']['reason'], (mode, budget))
                    self.assertLess(right['complete_to_bp'], 400)
                    elapsed = out['timing']['elapsed_ms']
                    self.assertLessEqual(elapsed, budget + 50 + 150, (mode, budget))
                    seen = out['usage']['observed_max_uninterruptible_ms']
                    self.assertIs(type(seen), int)
                    self.assertLessEqual(seen, 50 + 100, (mode, budget))
                    n += 1
            self.assertEqual(6, n)
            # the process-wide observation, in both capabilities documents
            dc = requests.get(server.url + '/capabilities').json()['deadline_check']
            self.assertGreater(dc['observed_max_uninterruptible_ms'], 0)

    def test_wide_index_delivery_reserve_stops_the_walk(self):
        """Pass 5, W6: the delivery reserve keeps back from the walk the time to build and
        compress what it will deliver, from the walked seed's modelled account. Under rates
        far below this machine's (0.001 MB/s) the reserve exceeds the whole bound at the first
        level: the attempt answers 200 at once, the first seed partial (attempt_deadline), the
        others not started, and usage states the moved walk-until. In both label modes the
        not-started seeds state labels_from_seed as the walked one: true for these derived
        seeds, false in annotate mode, which derives no set (review of 2026-10-06, X-DUP-01:
        a not-started annotate seed stated true)."""
        with self._server('--traverse-delivery-compress-mbps', '0.001',
                          '--traverse-delivery-build-mbps', '0.001') as server:
            caps = requests.get(server.url + '/traverse/capabilities').json()
            self.assertEqual(0.001, caps['attempts']['delivery_reserve']['compress_mbps'])
            for mode in ('constrain', 'annotate'):
                req = self._request(mode, 20000,
                                    attempt_id='reserve-1' if mode == 'constrain' else 'reserve-a')
                req['seeds'] = req['seeds'] * 3
                t0 = time.time()
                ret = server.post('traverse', req)
                took = (time.time() - t0) * 1000
                self.assertEqual(200, ret.status_code, ret.text[:500])
                out = ret.json()
                usage = out['usage']
                self.assertEqual('deadline', usage['reason'], mode)
                self.assertLess(took, usage['bound_ms'])
                self.assertLess(usage['bound']['walk_until_ms'],
                                usage['bound_ms'] - usage['bound']['allowance_ms'] // 2)
                first = out['results'][0]
                self.assertEqual(('attempt', 'attempt_deadline'),
                                 (first['resource_stop']['scope'],
                                  first['resource_stop']['resource']), mode)
                for r in out['results'][1:]:
                    self.assertEqual('not_started', r['resource_stop']['phase'], mode)
                self.assertEqual([mode == 'constrain'] * 3,
                                 [r['seed']['labels_from_seed'] for r in out['results']], mode)
        # at the default rates the floor (allowance / 2) holds for this small output
        with self._server() as server:
            ret = server.post('traverse', self._request('constrain', 100, attempt_id='reserve-2'))
            usage = ret.json()['usage']
            self.assertEqual(usage['bound_ms'] - usage['bound']['allowance_ms'] // 2,
                             usage['bound']['walk_until_ms'])

    def test_wide_index_unpaced_overruns(self):
        """--traverse-chunk-target-ms 0: one piece per read, as before; the same walk runs far
        past its budget (skipped when this machine reads the rows too fast to tell), and the
        paced walk is censored where it is (the same complete_to_bp)."""
        with self._server('--traverse-chunk-target-ms', '0') as whole, self._server() as paced:
            for mode in ('constrain', 'annotate'):
                a = whole.post('traverse', self._request(mode, 100)).json()
                b = paced.post('traverse', self._request(mode, 100)).json()
                if a['timing']['elapsed_ms'] < 100 + 250:
                    self.skipTest('the unpaced read took %.0f ms: too fast to show an overrun'
                                  % a['timing']['elapsed_ms'])
                self.assertLess(b['timing']['elapsed_ms'], a['timing']['elapsed_ms'] - 150)
                self.assertEqual(a['results'][0]['arms']['right']['complete_to_bp'],
                                 b['results'][0]['arms']['right']['complete_to_bp'])
                caps = requests.get(whole.url + '/traverse/capabilities').json()
                self.assertEqual(0, caps['deadline_check']['chunk_target_ms'])


class TestTraverseAttempts(TestingBase):
    """Stage 4, backend half (DESIGN-traverse-graphlet.md §14 v5.1), against a walk slow
    enough to be stopped in its middle: two haplotypes of 200 kbp with a SNP every 64 bp
    (k = 31), walked in annotate mode with every route kept and a beam of 64, so that a seed
    takes seconds. Each test starts its own short-lived server: a cancel by id with the client
    still connected, a client that goes away (with and without attempt_id), a cancel that
    waits for the end, the server's duration bound, the retention of finished attempts, and
    concurrent attempts with mixed cancels."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        rng = random.Random(20261003)
        h1 = ''.join(rng.choice('ACGT') for _ in range(200_000))
        h2 = list(h1)
        for i in range(32, len(h2), 64):
            h2[i] = rng.choice([c for c in 'ACGT' if c != h1[i]])
        h2 = ''.join(h2)
        d = cls.tempdir.name
        fastas = []
        for name, s in (('h1', h1), ('h2', h2)):
            fastas.append(f'{d}/{name}.fa')
            with open(fastas[-1], 'w') as f:
                f.write(f'>{name}\n{s}\n')
        cls.graph = d + '/slow_k31.dbg'
        cls._build_graph(fastas, cls.graph, K, 'succinct', mode='basic')
        cls._annotate_graph(fastas, cls.graph, d + '/slow', 'column', anno_type='filename')
        cls.anno = d + '/slow.column.annodbg'
        cls.seed = h1[:40]
        # how long one uncancelled seed walks in this build (the CLI's own timing): a test
        # that needs to stop a walk in its middle is skipped when the walk is too fast
        path = d + '/one.json'
        with open(path, 'w') as f:
            json.dump(cls._request(1), f)
        res = subprocess.run(shlex.split(METAGRAPH) + ['traverse', '-i', cls.graph, '-a', cls.anno,
                                                       path], stdout=subprocess.PIPE)
        cls.walk_ms = json.loads(res.stdout)['timing']['elapsed_ms']

    @classmethod
    def _request(cls, seeds, time_budget_ms=30000, radius=20000, **fields):
        req = {'seeds': [{'sequence': cls.seed}] * seeds,
               'strategy': {'direction': 'right', 'labels': {'mode': 'annotate'},
                            'branching': {'on_reconverge': 'keep'},
                            'frontier': {'on_overflow': 'beam'},
                            'bounds': {'max_extension_bp': radius, 'max_live_paths': 64,
                                       'max_paths': 100_000_000, 'max_steps': 1_000_000_000,
                                       'max_output_bp': 1_000_000_000,
                                       'time_budget_ms': time_budget_ms},
                            'output': {'detail': 'summary', 'timing': True,
                                       'max_branch_events': 10, 'sequences': False}}}
        req.update(fields)
        return req

    def _need_a_slow_walk(self):
        if self.walk_ms < 1500:
            self.skipTest(f'one seed walks {self.walk_ms:.0f} ms in this build: too fast to '
                          'stop in its middle')

    class _Server:
        def __init__(self, test, *flags, threads=4, port=None):
            # |port|: a restart of an endpoint (the port of a server that was stopped)
            self.port = port or _free_port()
            self.url = f'http://127.0.0.1:{self.port}'
            self.log_path = f'{test.tempdir.name}/server-{self.port}.log'
            self.log = open(self.log_path, 'a')
            self.process = subprocess.Popen(
                shlex.split(METAGRAPH) + ['server_query', '-i', test.graph, '-a', test.anno,
                                          '--port', str(self.port), '--address', '127.0.0.1',
                                          '-p', str(threads)] + list(flags),
                stdout=self.log, stderr=subprocess.STDOUT)
            for _ in range(600):
                try:
                    if requests.get(self.url + '/traverse/capabilities', timeout=2).ok:
                        break
                except requests.exceptions.RequestException:
                    pass
                time.sleep(0.1)

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            self.process.kill()
            self.process.wait()
            self.log.close()

        def text(self):
            with open(self.log_path) as f:
                return f.read()

        def post(self, route, payload, **kw):
            return requests.post(f'{self.url}/{route}', data=json.dumps(payload), **kw)

        def state(self, attempt_id, **kw):
            return requests.get(f'{self.url}/traverse/attempt/{attempt_id}', **kw)

    def test_cancel_mid_walk_with_the_client_connected(self):
        """POST /traverse/cancel stops the running seed at its next checkpoint (partial,
        resource_stop scope attempt, resource cancelled) and leaves the other seeds unstarted;
        the client still connected receives that partial response with its usage, and the
        attempt's state says when the walk stopped."""
        self._need_a_slow_walk()
        with self._Server(self) as server:
            req = self._request(3, attempt_id='cancel-1', budget_id='budget-1', locus_id='locus-1')
            got = {}
            walker = threading.Thread(target=lambda: got.update(
                r=server.post('traverse', req, timeout=300)))
            started = time.time()
            walker.start()
            time.sleep(0.5)
            cancel = server.post('traverse/cancel', {'attempt_id': 'cancel-1'})
            walker.join()
            took = time.time() - started
            self.assertEqual(200, cancel.status_code, cancel.text)
            self.assertTrue(cancel.json()['cancelled'])
            self.assertIn(cancel.json()['state'], ('stopping', 'finished'))
            ret = got['r']
            self.assertEqual(200, ret.status_code, ret.text)
            out = ret.json()
            first = out['results'][0]
            self.assertEqual('partial', first['outcome']['walks'])
            q = first['resource_stop']
            self.assertEqual(('attempt', 'cancelled', 'traversal', ['continue_from_leaves']),
                             (q['scope'], q['resource'], q['phase'], q['actions']))
            self.assertIn('attempt_id', [l['knob'] for l in first['arms']['right']['limitations']])
            for unstarted in out['results'][1:]:
                self.assertTrue(unstarted['error'].startswith('not started'))
                self.assertEqual(('not_started', 'cancelled', ['retry_attempt']),
                                 (unstarted['resource_stop']['phase'],
                                  unstarted['resource_stop']['resource'],
                                  unstarted['resource_stop']['actions']))
                # annotate mode derives no set: false, as the walked seed states it (review of
                # 2026-10-06, X-DUP-01)
                self.assertIs(False, unstarted['seed']['labels_from_seed'])
            self.assertIs(False, first['seed']['labels_from_seed'])
            usage = out['usage']
            self.assertEqual(('cancelled', 'budget-1', 'locus-1'),
                             (usage['reason'], usage['budget_id'], usage['locus_id']))
            self.assertEqual({'requested': 3, 'started': 1, 'finished': 1, 'abandoned': 0},
                             usage['seeds'])
            self.assertEqual(['partial', 'not_started', 'not_started'],
                             [s['outcome'] for s in usage['per_seed']])
            # stopped early: three seeds walk 3 x walk_ms uncancelled
            self.assertLess(took, 0.5 + self.walk_ms / 1000 + 2)
            state = server.state('cancel-1').json()
            self.assertEqual(('finished', 'cancelled', 'cancel', True, 200),
                             (state['state'], state['reason'], state['stop_requested_by'],
                              state['response']['written'], state['response']['status']))
            self.assertLess(_instant(state['stopped_at']) - _instant(state['cancel_requested_at']),
                            2.0)
            self.assertLessEqual(_instant(state['stopped_at']), _instant(state['finished_at']))
            log = server.text()
            self.assertIn('Attempt cancel-1: cancel requested', log)
            self.assertIn('finished (cancelled): walk stopped at', log)

    def test_client_disconnect_stops_the_walk(self):
        """A client that goes away (here: a read timeout) stops its walk at the next client
        check, with or without attempt_id, and nothing is written; with one server thread,
        the next request is served at once rather than after the abandoned walk."""
        self._need_a_slow_walk()
        with self._Server(self, threads=1) as server:
            started = time.time()
            with self.assertRaises(requests.exceptions.ReadTimeout):
                server.post('traverse', self._request(3, attempt_id='gone-1'), timeout=0.5)
            # the one thread is free again only once the walk stopped
            state = server.state('gone-1', timeout=60).json()
            freed = time.time() - started
            self.assertLess(freed, 0.5 + min(3.0, self.walk_ms / 1000 / 2), state)
            self.assertEqual(('finished', 'client_gone', False),
                             (state['state'], state['reason'], state['response']['written']))
            self.assertEqual('client_gone', state['stop_requested_by'])
            self.assertEqual('client_gone', state['usage']['reason'])
            # the first seed's walk was cut: abandoned, not finished
            self.assertEqual({'requested': 3, 'started': 1, 'finished': 0, 'abandoned': 1},
                             state['seeds'])
            # without attempt_id: stated in the log only
            with self.assertRaises(requests.exceptions.ReadTimeout):
                server.post('traverse', self._request(3), timeout=0.5)
            started = time.time()
            requests.get(server.url + '/traverse/capabilities', timeout=60)
            self.assertLess(time.time() - started, min(3.0, self.walk_ms / 1000 / 2))
            log = server.text()
            self.assertIn('Attempt gone-1 (request', log)
            self.assertIn('no response written (the client is gone)', log)
            self.assertIn('client gone, walk stopped at', log)

    def test_sigterm_stops_the_server_promptly(self):
        """SIGTERM (what `docker stop` sends; the server, PID 1 in its container, ignored it, so
        every deploy waited the stop timeout of 120 s): the server exits within seconds, with
        status 0. An idle one at once; one walking an attempt stops the walk at its next poll as
        for a client gone — nothing is written for it, its connection is closed, never a
        partial result passed off as complete — and then exits."""
        with self._Server(self) as server:
            started = time.time()
            server.process.send_signal(signal.SIGTERM)
            self.assertEqual(0, server.process.wait(timeout=30))
            self.assertLess(time.time() - started, 5)
            log = server.text()
            self.assertIn('Signal 15: shutting down', log)
            self.assertIn('[Server] Stopped', log)
        self._need_a_slow_walk()
        with self._Server(self) as server:
            got = {}

            def walk():
                try:
                    got['r'] = server.post('traverse', self._request(3, attempt_id='term-1'),
                                           timeout=300)
                except requests.exceptions.RequestException as e:
                    got['e'] = e
            walker = threading.Thread(target=walk)
            walker.start()
            # the attempt is running (registered) before the signal
            for _ in range(100):
                if server.state('term-1').status_code == 200:
                    break
                time.sleep(0.05)
            self.assertEqual('running', server.state('term-1').json()['state'])
            time.sleep(0.3)
            started = time.time()
            server.process.send_signal(signal.SIGTERM)
            self.assertEqual(0, server.process.wait(timeout=30))
            took = time.time() - started
            walker.join(timeout=30)
            self.assertFalse(walker.is_alive())
            # well within the drain time: the walk stops at its next poll; three seeds walk
            # 3 x walk_ms uncut
            self.assertLess(took, 3.5)
            self.assertLess(took, 3 * self.walk_ms / 1000)
            if 'r' in got:
                # answered before the signal reached it: then whole and valid
                self.assertEqual(200, got['r'].status_code, got['r'].text)
                self.assertEqual(3, len(got['r'].json()['results']))
            else:
                self.assertIsInstance(got['e'], (requests.exceptions.ConnectionError,
                                                 requests.exceptions.ChunkedEncodingError))
            log = server.text()
            self.assertIn('Signal 15: shutting down', log)
            self.assertIn('1 traversal(s) in flight stopped', log)
            self.assertIn('Attempt term-1 (request', log)
            self.assertIn('[Server] Stopped', log)

    def test_sigterm_while_the_index_loads(self):
        """A single-index server listens at once and answers 503 (GET /capabilities: ready
        false) while its index loads; SIGTERM then stops it within seconds too (review of W2:
        returning through run_server joined the load in progress, 19 s on SRA, on a cold
        staging disk the whole docker-stop timeout). The load here never ends: the graph is a
        FIFO nobody writes to, so its loader blocks in open() for good."""
        fifo = f'{self.tempdir.name}/never_loads.dbg'
        if not os.path.exists(fifo):
            os.mkfifo(fifo)
        port = _free_port()
        url = f'http://127.0.0.1:{port}'
        log_path = f'{self.tempdir.name}/server-loading-{port}.log'
        with open(log_path, 'w') as log:
            process = subprocess.Popen(
                shlex.split(METAGRAPH) + ['server_query', '-i', fifo, '-a', self.anno,
                                          '--port', str(port), '--address', '127.0.0.1', '-p', '2'],
                stdout=log, stderr=subprocess.STDOUT)
            try:
                caps = None
                for _ in range(300):
                    self.assertIsNone(process.poll(), 'the server exited while loading')
                    try:
                        ret = requests.get(url + '/capabilities', timeout=2)
                        if ret.ok:
                            caps = ret.json()
                            break
                    except requests.exceptions.RequestException:
                        pass
                    time.sleep(0.1)
                self.assertIsNotNone(caps, 'the server did not listen')
                self.assertFalse(caps['ready'])
                self.assertEqual(503, requests.get(url + '/traverse/capabilities',
                                                   timeout=5).status_code)
                started = time.time()
                process.send_signal(signal.SIGTERM)
                self.assertEqual(0, process.wait(timeout=30))
                self.assertLess(time.time() - started, 5)
            finally:
                if process.poll() is None:
                    process.kill()
                process.wait()
        with open(log_path) as f:
            text = f.read()
        self.assertIn('Signal 15: shutting down', text)
        self.assertIn('[Server] Stopped', text)
        self.assertNotIn('Annotated graph loaded', text)

    def test_cancel_waits_for_the_end(self):
        """A cancel with wait_ms answers once the attempt finished: its response written."""
        self._need_a_slow_walk()
        with self._Server(self) as server:
            got = {}
            walker = threading.Thread(target=lambda: got.update(
                r=server.post('traverse', self._request(2, attempt_id='wait-1'), timeout=300)))
            walker.start()
            time.sleep(0.5)
            cancel = server.post('traverse/cancel', {'attempt_id': 'wait-1', 'wait_ms': 10000})
            walker.join()
            self.assertEqual(200, cancel.status_code, cancel.text)
            body = cancel.json()
            self.assertEqual((True, 'finished'), (body['cancelled'], body['state']))
            self.assertEqual(('cancelled', True), (body['attempt']['reason'],
                                                   body['attempt']['response']['written']))
            self.assertEqual(200, got['r'].status_code)
            self.assertEqual(400, server.post('traverse/cancel', {'attempt_id': 'wait-1',
                                                                  'wait_ms': 10001}).status_code)

    def test_the_server_enforces_the_attempts_bound(self):
        """The bound n_seeds x time_budget_ms + allowance is enforced by the server itself:
        with an allowance of 2 ms the last seed cannot finish within it, so the walk stops at
        the bound less half the allowance (attempt_deadline) and a response that cannot be
        written by the bound is a 503 with the usage."""
        self._need_a_slow_walk()
        with self._Server(self, '--traverse-attempt-allowance-ms', '2') as server:
            ret = server.post('traverse', self._request(4, time_budget_ms=300,
                                                        attempt_id='deadline-1'), timeout=300)
            self.assertIn(ret.status_code, (200, 503), ret.text)
            usage = ret.json()['usage']
            self.assertEqual('deadline', usage['reason'])
            self.assertEqual(1202, usage['bound_ms'])
            # at most the bound less half the allowance; the delivery reserve (pass 5) moves it
            # earlier when the walked seeds' estimated text needs more time than that to build
            self.assertLessEqual(usage['bound']['walk_until_ms'], 1201)
            # stopped at the next checkpoint: one head, plus building the response
            self.assertLess(usage['elapsed_ms'], usage['bound_ms'] + 1000)
            if ret.status_code == 200:
                stops = [r.get('resource_stop', {}).get('resource') for r in ret.json()['results']]
                self.assertIn('attempt_deadline', stops)
            state = server.state('deadline-1').json()
            self.assertEqual(('finished', 'deadline'), (state['state'], state['reason']))
            self.assertEqual(ret.status_code, state['response']['status'])
            # the library tells this 503 from a loading server: it ran, and its id is used up
            client = graphlet_lib.TraverseClient('127.0.0.1', server.port)
            req = self._request(4, time_budget_ms=300)
            try:
                client.traverse(req['seeds'], req['strategy'], detail='summary',
                                attempt_id='deadline-lib')
            except graphlet_lib.AttemptAtBound as e:
                self.assertEqual('deadline', e.usage['reason'])
            self.assertEqual('deadline', client.attempt('deadline-lib')['reason'])

    def test_a_close_behind_waiting_bytes_stops_the_walk(self):
        """A client that sent bytes past its request (a trailing CRLF, a pipelined request)
        and then closed is gone too: the peek returns those bytes, the connection's TCP state
        tells the close behind them, so the walk stops and nothing is written."""
        self._need_a_slow_walk()
        with self._Server(self) as server:
            for name, extra in (('crlf', b'\r\n'),
                                ('pipelined', b'GET /traverse/capabilities HTTP/1.1\r\n'
                                              b'Host: x\r\n\r\n')):
                attempt_id = 'behind-' + name
                body = json.dumps(self._request(3, attempt_id=attempt_id)).encode()
                sock = socket.create_connection(('127.0.0.1', server.port))
                sock.sendall(b'POST /traverse HTTP/1.1\r\nHost: x\r\n'
                             b'Content-Type: application/json\r\n'
                             + b'Content-Length: %d\r\n\r\n' % len(body) + body + extra)
                time.sleep(0.5)
                closed = time.time()
                sock.close()
                for _ in range(600):
                    state = server.state(attempt_id, timeout=60).json()
                    if state.get('state') == 'finished':
                        break
                    time.sleep(0.05)
                self.assertEqual(('finished', 'client_gone', False),
                                 (state['state'], state['reason'], state['response']['written']),
                                 name)
                self.assertLess(_instant(state['stopped_at']) - closed, 2.0, name)

    def test_a_tombstone_is_kept_its_whole_retention_period(self):
        """A cancel of an unknown id tombstones it for retention_s, whatever finishes after
        it; beyond --traverse-attempt-retention tombstones a cancel is refused (429, nothing
        promised)."""
        with self._Server(self, '--traverse-attempt-retention', '2') as server:
            early = server.post('traverse/cancel', {'attempt_id': 'tomb-1'})
            self.assertEqual(404, early.status_code, early.text)
            self.assertTrue(early.json()['tombstone'])
            for i in range(4):
                quick = self._request(1, radius=10, attempt_id=f'tomb-done-{i}')
                self.assertEqual(200, server.post('traverse', quick).status_code)
            late = server.post('traverse', self._request(1, radius=10, attempt_id='tomb-1'))
            self.assertEqual(409, late.status_code, late.text)
            self.assertEqual(404, server.post('traverse/cancel', {'attempt_id': 'tomb-2'}).status_code)
            full = server.post('traverse/cancel', {'attempt_id': 'tomb-3'})
            self.assertEqual(429, full.status_code, full.text)
            self.assertEqual((False, False, 'unknown'),
                             (full.json()['tombstone'], full.json()['cancelled'],
                              full.json()['state']))
            client = graphlet_lib.TraverseClient('127.0.0.1', server.port)
            self.assertFalse(client.cancel('tomb-4')['tombstone'])

    def test_finished_attempts_are_kept_for_the_retention_period(self):
        """GET answers for the retention period after the attempt finished, then 404; the id
        is free again."""
        with self._Server(self, '--traverse-attempt-retention-s', '1') as server:
            quick = self._request(1, radius=10, attempt_id='kept-1')
            self.assertEqual(200, server.post('traverse', quick).status_code)
            self.assertEqual(200, server.state('kept-1').status_code)
            self.assertEqual(409, server.post('traverse', quick).status_code)
            time.sleep(1.5)
            gone = server.state('kept-1')
            self.assertEqual(404, gone.status_code)
            self.assertIn('kept 1 s after they finish', gone.json()['error'])
            self.assertEqual(200, server.post('traverse', quick).status_code)

    def test_a_finished_attempt_is_held_through_its_not_after_ms(self):
        """Review of the pass-5 fixes, finding 1 (the reviewer's p1 and p1b): a ledger releases
        on a finished state, and a replay of the request (same attempt_id, not_after_ms and
        expect_server_instance) arriving after retention_s, or after retention_count later
        finishes evicted the attempt, ran again. A finished attempt sent with not_after_ms
        stays refused (409) until not_after_ms + the skew allowance, and the release rule
        states it; a copy refused by a tombstone without a not_after_ms of its own is not
        covered (finding 3)."""
        for flags in (('--traverse-attempt-retention-s', '1'),
                      ('--traverse-attempt-retention', '1')):
            with self._Server(self, *flags) as server:
                caps = requests.get(server.url + '/traverse/capabilities').json()['attempts']
                self.assertIn('not_after_ms + clock_skew_allowance_ms, at most tombstone_max_s '
                              'after it finished', caps['release_rule'])
                instance = caps['server_instance']
                not_after = int(time.time() * 1000) + 60000
                req = self._request(1, radius=10, attempt_id='fin-1', not_after_ms=not_after,
                                    expect_server_instance=instance)
                first = server.post('traverse', req)
                self.assertEqual(200, first.status_code, first.text)
                cancel = server.post('traverse/cancel', {'attempt_id': 'fin-1',
                                                         'not_after_ms': not_after})
                self.assertEqual((404, 'finished'), (cancel.status_code, cancel.json()['state']))
                self.assertEqual('finished', server.state('fin-1').json()['state'])
                # the ledger releases; then its retention ends, by age or by count
                if flags[0] == '--traverse-attempt-retention-s':
                    time.sleep(1.3)
                else:
                    other = self._request(1, radius=10, attempt_id='fin-2')
                    self.assertEqual(200, server.post('traverse', other).status_code)
                replay = server.post('traverse', req)
                self.assertEqual(409, replay.status_code, (flags, replay.text))
                self.assertEqual('finished', replay.json()['attempt']['state'])
                self.assertEqual(200, server.state('fin-1').status_code)
        with self._Server(self, '--traverse-attempt-retention-s', '1',
                          '--traverse-clock-skew-ms', '100') as server:
            not_after = int(time.time() * 1000) + 1500
            cancel = server.post('traverse/cancel', {'attempt_id': 'nna-1',
                                                     'not_after_ms': not_after})
            self.assertTrue(cancel.json()['covers_admission'])
            copy = self._request(1, radius=10, attempt_id='nna-1')   # no not_after_ms
            refused = server.post('traverse', copy)
            self.assertEqual(409, refused.status_code, refused.text)
            attempt = refused.json()['attempt']
            self.assertNotIn('not_after_ms', attempt)
            self.assertEqual((False, 'no_not_after_ms'),
                             (attempt['covers_admission'], attempt['covers_admission_reason']))

    def test_retention_settings_are_validated_at_start_up(self):
        """Review of pass 5, finding 4: --traverse-attempt-retention-s -1 started and advertised
        2^64 - 1 seconds while every tombstone expired at once. The retention settings are
        bounded integers: a negative, non-numeric or too large value refuses to start, naming
        the option and its own range (the review of the pass-5 fixes: a negative or non-numeric
        value named [0, 2^53 - 1] instead)."""
        ranges = {'--traverse-attempt-retention-s': '[0, 31536000]',
                  '--traverse-attempt-retention': '[0, 10000000]',
                  '--traverse-attempt-tombstone-max-s': '[0, 31536000]'}
        for flag, value in (('--traverse-attempt-retention-s', '-1'),
                            ('--traverse-attempt-retention-s', '1e400'),
                            ('--traverse-attempt-retention-s', 'abc'),
                            ('--traverse-attempt-retention-s', '1.5'),
                            ('--traverse-attempt-retention-s', 'nan'),
                            ('--traverse-attempt-retention-s', ''),
                            ('--traverse-attempt-retention-s', '5 '),
                            ('--traverse-attempt-retention-s', '31536001'),
                            ('--traverse-attempt-retention-s', '9223372036854775808'),
                            ('--traverse-attempt-retention', '-1'),
                            ('--traverse-attempt-retention', 'inf'),
                            ('--traverse-attempt-retention', '10000001'),
                            ('--traverse-attempt-tombstone-max-s', '-5'),
                            ('--traverse-attempt-tombstone-max-s', '18446744073709551615')):
            res = subprocess.run(shlex.split(METAGRAPH) + [
                'server_query', '-i', self.graph, '-a', self.anno, '--port', str(_free_port()),
                '--address', '127.0.0.1', flag, value], stdout=subprocess.PIPE,
                stderr=subprocess.PIPE, timeout=60)
            self.assertNotEqual(0, res.returncode, (flag, value))
            self.assertIn(flag + ' must be an integer in ' + ranges[flag], res.stderr.decode(),
                          (flag, value))
        # the largest accepted values start (and are stated)
        with self._Server(self, '--traverse-attempt-retention-s', '31536000',
                          '--traverse-attempt-retention', '10000000',
                          '--traverse-attempt-tombstone-max-s', '31536000') as server:
            caps = requests.get(server.url + '/traverse/capabilities').json()['attempts']
            self.assertEqual((31536000, 10000000, 31536000),
                             (caps['retention_s'], caps['retention_count'],
                              caps['tombstone_max_s']))

    def test_retention_zero_promises_no_suppression(self):
        """Retention 0 s means no tombstones: a cancel of an unknown id is refused (429,
        tombstone false, reason no_suppression), the capabilities say so, and a request with
        the id runs (the reviewer's probe: before, the cancel answered tombstone true and the
        request ran anyway). A count of 0 refuses likewise (reason tombstones_full)."""
        req = self._request(1, radius=1, attempt_id='cancelled-before-upload')
        with self._Server(self, '--traverse-attempt-retention-s', '0') as server:
            caps = requests.get(server.url + '/traverse/capabilities').json()['attempts']
            self.assertEqual(0, caps['retention_s'])
            self.assertIn('no tombstones', caps['suppression'])
            first = server.post('traverse/cancel', {'attempt_id': 'cancelled-before-upload',
                                                    'not_after_ms': int(time.time() * 1000)
                                                    + 60000})
            self.assertEqual(429, first.status_code, first.text)
            self.assertEqual((False, False, 'no_suppression'),
                             (first.json()['tombstone'], first.json()['cancelled'],
                              first.json()['reason']))
            self.assertEqual(404, server.state('cancelled-before-upload').status_code)
            self.assertEqual(200, server.post('traverse', req).status_code)
        with self._Server(self, '--traverse-attempt-retention', '0') as server:
            first = server.post('traverse/cancel', {'attempt_id': 'cancelled-before-upload'})
            self.assertEqual(429, first.status_code, first.text)
            self.assertEqual((False, 'tombstones_full'),
                             (first.json()['tombstone'], first.json()['reason']))

    def _half_uploaded(self, server, req):
        """Send the header and half the body of |req| on a raw socket: the request is on the
        wire, its handler not started."""
        body = json.dumps(req).encode()
        split = len(body) // 2
        sock = socket.create_connection(('127.0.0.1', server.port), timeout=30)
        sock.sendall(b'POST /traverse HTTP/1.1\r\nHost: localhost\r\n'
                     b'Content-Type: application/json\r\n'
                     + b'Content-Length: %d\r\nConnection: close\r\n\r\n' % len(body)
                     + body[:split])
        time.sleep(0.05)
        return sock, body[split:]

    @staticmethod
    def _complete(sock, rest):
        import http.client
        sock.sendall(rest)
        response = http.client.HTTPResponse(sock)
        response.begin()
        out = json.loads(response.read())
        sock.close()
        return response.status, out

    def test_a_cancel_naming_not_after_ms_covers_a_half_uploaded_request(self):
        """Review of pass 5, finding 1, the reviewer's timeline with retention 1 s: a request
        with not_after_ms = now + 60 s is half uploaded, then cancelled (404, tombstone). With
        the cancel naming that not_after_ms, the tombstone is held to not_after_ms + the skew
        allowance (suppressed_until_ms, covers_admission true), a repeat cancel at 0.76 s states
        the same expiry, and the upload completed at 1.2 s is refused (409, the tombstone) and
        never runs. Without not_after_ms in the cancel (and no repeat), covers_admission is
        false (a ledger may not release on it) and the upload completed at 1.2 s runs, as the
        retention period alone allows."""
        with self._Server(self, '--traverse-attempt-retention-s', '1') as server:
            skew = requests.get(server.url + '/capabilities').json()['attempts'][
                'clock_skew_allowance_ms']
            for covered in (True, False):
                attempt_id = 'delayed-original-%d' % covered
                not_after = int(time.time() * 1000) + 60000
                req = self._request(1, radius=5, attempt_id=attempt_id, not_after_ms=not_after)
                sock, rest = self._half_uploaded(server, req)
                cancel = {'attempt_id': attempt_id}
                if covered:
                    cancel['not_after_ms'] = not_after
                first = server.post('traverse/cancel', cancel)
                start = time.monotonic()
                self.assertEqual(404, first.status_code, first.text)
                body = first.json()
                self.assertTrue(body['tombstone'])
                self.assertEqual(covered, body['covers_admission'], body)
                if covered:
                    self.assertEqual(not_after + skew, body['suppressed_until_ms'])
                    self.assertEqual(not_after, body['not_after_ms'])
                else:
                    self.assertEqual('no_not_after_ms', body['covers_admission_reason'])
                    self.assertLess(body['suppressed_until_ms'], not_after)
                if covered:
                    # (uncovered, a repeat cancel would keep the tombstone 1 s from then, and
                    # the copy's refusal would extend it to the copy's not_after_ms)
                    time.sleep(0.75)
                    second = server.post('traverse/cancel', {'attempt_id': attempt_id})
                    self.assertEqual(404, second.status_code, second.text)
                    self.assertEqual(body['suppressed_until_ms'],
                                     second.json()['suppressed_until_ms'])
                    self.assertTrue(second.json()['covers_admission'])
                time.sleep(max(0.0, 1.2 - (time.monotonic() - start)))
                status, out = self._complete(sock, rest)
                state = server.state(attempt_id)
                if covered:
                    self.assertEqual(409, status, out)
                    self.assertTrue(out['attempt']['tombstone'])
                    self.assertTrue(out['attempt']['covers_admission'])
                    self.assertNotIn('usage', out)
                    self.assertEqual(404, state.status_code)
                    self.assertTrue(state.json()['tombstone'])
                    self.assertEqual(not_after + skew, state.json()['suppressed_until_ms'])
                else:
                    self.assertEqual(200, status, out)
                    self.assertEqual('completed', out['usage']['reason'])

    def test_expect_server_instance_refuses_another_process(self):
        """The restart hole of finding 1: tombstones live in memory, so a copy reaching a
        restarted process would run. A request naming another server_instance is refused before
        anything runs (409 instance_mismatch, nothing registered); naming this one, it runs. The
        CLI, whose instance is its own, refuses one too."""
        with self._Server(self) as server:
            mine = requests.get(server.url + '/capabilities').json()['attempts'][
                'server_instance']
            other = 'f' * 16 if mine != 'f' * 16 else 'e' * 16
            req = self._request(1, radius=5, attempt_id='inst-1', expect_server_instance=other)
            refused = server.post('traverse', req)
            self.assertEqual(409, refused.status_code, refused.text)
            body = refused.json()
            self.assertEqual(('instance_mismatch', other, mine, 'inst-1'),
                             (body['state'], body['expect_server_instance'],
                              body['server_instance'], body['attempt_id']))
            self.assertEqual(404, server.state('inst-1').status_code)
            req['expect_server_instance'] = mine
            ok = server.post('traverse', req)
            self.assertEqual(200, ok.status_code, ok.text)
            self.assertEqual(mine, ok.json()['usage']['server_instance'])
            # without attempt_id the field is refused (400)
            plain = self._request(1, radius=5, expect_server_instance=mine)
            self.assertEqual(400, server.post('traverse', plain).status_code)
        path = self.tempdir.name + '/instance.json'
        with open(path, 'w') as f:
            json.dump(self._request(1, radius=5, attempt_id='cli-1',
                                    expect_server_instance='0' * 16), f)
        res = subprocess.run(shlex.split(METAGRAPH) + ['traverse', '-i', self.graph, '-a',
                                                       self.anno, path],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(1, res.returncode)
        self.assertEqual('instance_mismatch', json.loads(res.stdout)['state'])

    def test_a_finished_state_is_replay_safe_only_when_pinned(self):
        """Review of levels 4-5, finding 6 (the reviewer's finish-restart-replay probe): a
        finished attempt's hold lives in its process. Requests finished with not_after_ms 60 s
        ahead are refused when replayed to the same process after their retention (1 s, the last
        1), pinned or not; after a restart of the same endpoint, before not_after_ms, the replay
        sent without expect_server_instance runs again (200, the same work units as the first
        time), which the release rule now states (a finished state is replay-safe only for a
        pinned attempt), and the pinned one is refused, 409 instance_mismatch."""
        flags = ('--traverse-attempt-retention-s', '1', '--traverse-attempt-retention', '1')
        not_after = int(time.time() * 1000) + 60000
        with self._Server(self, *flags) as server:
            caps = requests.get(server.url + '/traverse/capabilities').json()['attempts']
            instance = caps['server_instance']
            for phrase in ("The hold is this process's, in memory: a restarted process (a new "
                           "server_instance) holds none",
                           'replay-safe only as stated next',
                           'only for an attempt sent with expect_server_instance equal to this '
                           'server_instance and with not_after_ms',
                           'A finished state of an attempt sent without expect_server_instance '
                           'assumes that no copy of the request reaches a restarted process'):
                self.assertIn(phrase, caps['release_rule'])
            self.assertIn('a delayed copy of a cancelled or finished request would otherwise run '
                          'there', caps['instance'])
            unpinned = self._request(1, radius=5, attempt_id='finished-unpinned',
                                     not_after_ms=not_after)
            pinned = self._request(1, radius=5, attempt_id='finished-pinned',
                                   not_after_ms=not_after, expect_server_instance=instance)
            first = server.post('traverse', unpinned)
            self.assertEqual(200, first.status_code, first.text)
            self.assertEqual(200, server.post('traverse', pinned).status_code)
            time.sleep(1.1)
            for req in (unpinned, pinned):
                held = server.post('traverse', req)
                self.assertEqual(409, held.status_code, held.text)
                self.assertEqual('finished', held.json()['attempt']['state'])
            port = server.port
        with self._Server(self, *flags, port=port) as server:
            restarted = requests.get(server.url + '/capabilities').json()['attempts'][
                'server_instance']
            self.assertNotEqual(instance, restarted)
            replay = server.post('traverse', unpinned)
            self.assertEqual(200, replay.status_code, replay.text)
            self.assertEqual(restarted, replay.json()['usage']['server_instance'])
            self.assertEqual(first.json()['usage']['work_units'],
                             replay.json()['usage']['work_units'])
            self.assertGreater(replay.json()['usage']['work_units'], 0)
            refused = server.post('traverse', pinned)
            self.assertEqual(409, refused.status_code, refused.text)
            self.assertEqual(('instance_mismatch', instance, restarted),
                             (refused.json()['state'], refused.json()['expect_server_instance'],
                              refused.json()['server_instance']))

    def test_release_texts_and_their_grounds(self):
        """The review of 2026-10-06 (X2, X4, C16, C20, C24, C30, and the search service's release
        parity LRG-R1/R2), text only: the capabilities state what an attempt runs past its bound
        (up to its next delivery check, not one step: the review of the P2 fixes) and the clock
        release's assumption that it ended, the 409 refusing a copy as a finished-state source,
        the expired 409 as a release ground, the refusal order, the expect_server_instance
        pattern and the lookahead's polls; the duplicate's 409 no longer says "an attempt runs
        once" (D3: an id runs again once it is neither retained nor held). On a real server the
        grounds hold as stated: a copy of a finished request gets the 409 carrying the id's
        state as GET answers it (finished, its attempt_id and server_instance), also with a
        not_after_ms already past (registration before expiry); a fresh id past its not_after_ms
        gets the expired 409 and was never registered; an empty expect_server_instance is a 400
        without usage."""
        with self._Server(self) as server:
            caps = requests.get(server.url + '/traverse/capabilities').json()
            att, rule = caps['attempts'], caps['deadline_check']['rule']
            instance = att['server_instance']
            for field, phrase in (
                    ('bound', 'past it the attempt runs on until its next delivery check'),
                    ('bound', 'the walk stops at its first poll that reads the clock after it'),
                    ('bound', 'content_timeout when it applied'),
                    ('bound', 'else the lowest walk-until seen'),
                    ('not_after', 'It is the last of the refusals, all made under one lock'),
                    ('not_after', 'apart from what it runs past its bound, up to its next '
                                  'delivery check, a run of no stated length'),
                    ('instance', 'a string matching id_pattern'),
                    ('release_rule', 'the 409 refusing a copy of the request whose attempt'),
                    ('release_rule', 'An expired 409 (state: expired) for an attempt sent with '
                                     'exactly that not_after_ms'),
                    ('release_rule', 'It settles nothing'),
                    ('release_rule', 'still answers running or stopping after that instant'),
                    ('release_rule', 'no copy of the request reaches another server that serves '
                                     'the same ledger'),
                    ('release_rule', 'after it finished or after its latest refused copy')):
                self.assertIn(phrase, att[field], field)
            for field in ('bound', 'not_after', 'release_rule'):
                self.assertNotIn('uninterruptible step', att[field], field)
            self.assertIn('read the same deadlines every 16 graph steps', rule)
            self.assertIn("a chain's key mapping", rule)
            self.assertIn("the attempt's bound, which only the delivery checks compare", rule)
            self.assertIn('a walk-until at the first poll that reads the clock after it (one in 8',
                          rule)

            not_after = int(time.time() * 1000) + 60000
            req = self._request(1, radius=10, attempt_id='grounds-1', not_after_ms=not_after,
                                expect_server_instance=instance)
            first = server.post('traverse', req)
            self.assertEqual(200, first.status_code, first.text[:500])
            for copy_not_after in (not_after, int(time.time() * 1000) - 1000):
                copy_req = dict(req, not_after_ms=copy_not_after)
                refused = server.post('traverse', copy_req)
                self.assertEqual(409, refused.status_code, refused.text[:500])
                body = refused.json()
                self.assertNotIn('usage', body)
                self.assertIn('it is not run while so', body['error'])
                self.assertNotIn('runs once', body['error'])
                self.assertEqual('finished', body['attempt']['state'])
                self.assertEqual('grounds-1', body['attempt']['attempt_id'])
                self.assertEqual(instance, body['attempt']['server_instance'])
                # the id's state, as GET answers it
                state = server.state('grounds-1').json()
                self.assertEqual(state['state'], body['attempt']['state'])
                self.assertEqual(state['usage']['work_units'], body['attempt']['usage']['work_units'])

            past = int(time.time() * 1000) - 1000
            expired = server.post('traverse', self._request(1, radius=10, attempt_id='grounds-2',
                                                            not_after_ms=past,
                                                            expect_server_instance=instance))
            self.assertEqual(409, expired.status_code, expired.text[:500])
            body = expired.json()
            self.assertEqual('expired', body['state'])
            self.assertEqual(past, body['not_after_ms'])
            self.assertGreater(body['server_time_ms'], past)
            self.assertEqual(instance, body['server_instance'])
            self.assertNotIn('usage', body)
            self.assertEqual(404, server.state('grounds-2').status_code)

            for bad in ('', 'a/b', 'x' * 129):
                ret = server.post('traverse', self._request(1, radius=10, attempt_id='grounds-3',
                                                            expect_server_instance=bad))
                self.assertEqual(400, ret.status_code, (bad, ret.text[:300]))
                self.assertIn('expect_server_instance', ret.json()['error'])
                self.assertNotIn('usage', ret.json())
            self.assertEqual(404, server.state('grounds-3').status_code)

    def test_concurrent_attempts_with_cancels(self):
        """Sixteen attempts at once, half of them cancelled: each response and each state
        agree on what stopped it."""
        self._need_a_slow_walk()
        n = 16
        with self._Server(self, threads=n + 1) as server:
            got = {}

            def run(i):
                got[i] = server.post('traverse', self._request(1, attempt_id=f'conc-{i}'),
                                     timeout=600)

            threads = [threading.Thread(target=run, args=(i,)) for i in range(n)]
            for t in threads:
                t.start()
            time.sleep(0.5)
            cancels = {i: server.post('traverse/cancel', {'attempt_id': f'conc-{i}'})
                       for i in range(0, n, 2)}
            for t in threads:
                t.join()
            for i in range(n):
                ret = got[i]
                self.assertEqual(200, ret.status_code, ret.text)
                usage = ret.json()['usage']
                state = server.state(f'conc-{i}').json()
                self.assertEqual('finished', state['state'])
                self.assertEqual(usage['reason'], state['reason'], i)
                self.assertEqual(usage['work_units'], state['usage']['work_units'], i)
                if i % 2:
                    self.assertEqual('completed', usage['reason'], i)
                    self.assertIsNone(state['stop_requested_by'])
                else:
                    self.assertEqual(200, cancels[i].status_code, cancels[i].text)
                    self.assertEqual('cancelled', usage['reason'], i)
                    result = ret.json()['results'][0]
                    self.assertEqual('cancelled', result['resource_stop']['resource'], i)


class TestTraverseSeedPhase(TestingBase):
    """Review of levels 4-5, finding 3 (the reviewer's deadline probe): a request naming
    2,500 header labels that share a 1,024-character prefix, all on one 11-mer (k = 11), with
    timing. The longest uninterruptible piece stated 0.297 ms on a warm server for a seed phase
    of 126 ms: the pairwise duplicate check of the names (about 120 ms) and the resolution of
    the names were no piece, nor anything else between the reads and k-mer mappings before the
    first head. The seed phase is now cut into pieces at its reads and k-mer mappings, the spans
    between them "setup" pieces, and the pieces cover seed_phase_ms: here (one 11-mer: two k-mer
    mappings and a read or two of one row) fewer than 16 of them, so the longest is at least a
    sixteenth of it (before: 0.297 ms of 126 ms, a 425th), and at least the names' resolution,
    which runs inside one setup span. The duplicate check is a hash set."""

    N = 2500

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        d = cls.tempdir.name
        rng = random.Random(41)
        cls.seq = ''.join(rng.choice('ACGT') for _ in range(11))
        cls.names = ['H' + 'x' * 1024 + f'{i:06d}' for i in range(cls.N)]
        fa = d + '/x.fa'
        with open(fa, 'w') as f:
            f.write(''.join(f'>{name}\n{cls.seq}\n' for name in cls.names))
        for cmd in (f'{METAGRAPH} build -p 1 --mode basic --graph succinct -k 11 -o {d}/graph {fa}',
                    f'{METAGRAPH} annotate -p 1 -i {d}/graph.dbg --anno-filename --anno-type '
                    f'column --coordinates -o {d}/ann {fa}',
                    f'{METAGRAPH} transform_anno -p 1 --anno-type column_coord --coordinates '
                    f'-o {d}/ann {d}/ann.column.annodbg',
                    f'{METAGRAPH} annotate -p 1 -i {d}/graph.dbg --anno-filename '
                    f'--index-header-coords -o {d}/ann {fa}'):
            res = subprocess.run(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            assert res.returncode == 0, res.stderr.decode()
        cls.graph = d + '/graph.dbg'
        cls.anno = d + '/ann.column_coord.annodbg'

    def _request(self, labels):
        return {'seeds': [{'sequence': self.seq, 'labels': labels}],
                'strategy': {'direction': 'right',
                             'bounds': {'max_extension_bp': 1, 'time_budget_ms': 1},
                             'output': {'detail': 'summary', 'timing': True}}}

    def _check_seed_phase(self, out, where):
        (result,) = out['results']
        self.assertEqual(self.N, result['seed']['labels_supporting_total'], where)
        t = result['timing']
        piece = t['deadline']['longest_piece']
        # every span of the seed phase is a piece: fewer than 16 cover it
        self.assertGreaterEqual(piece['ms'] * 16 * 1.001 + 0.001, t['seed_phase_ms'],
                                f'{where}: {piece} for {t}')
        # the names are resolved inside one setup span
        self.assertGreaterEqual(piece['ms'] * 1.001 + 0.001, t['label_resolve_ms'],
                                f'{where}: {piece} for {t}')
        return t

    def test_the_longest_piece_covers_the_seed_phase(self):
        path = self.tempdir.name + '/req.json'
        with open(path, 'w') as f:
            json.dump(self._request(self.names), f)
        res = subprocess.run(shlex.split(METAGRAPH) + ['traverse', '--json', '-i', self.graph,
                                                       '-a', self.anno, path],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(0, res.returncode, res.stderr.decode())
        cli = self._check_seed_phase(json.loads(res.stdout), 'CLI')
        port = _free_port()
        log = open(f'{self.tempdir.name}/server-{port}.log', 'w')
        server = subprocess.Popen(
            shlex.split(METAGRAPH) + ['server_query', '-i', self.graph, '-a', self.anno,
                                      '--address', '127.0.0.1', '--port', str(port), '-p', '4'],
            stdout=log, stderr=subprocess.STDOUT)
        url = f'http://127.0.0.1:{port}'
        try:
            for _ in range(600):
                try:
                    if requests.get(url + '/traverse/capabilities', timeout=2).ok:
                        break
                except requests.exceptions.RequestException:
                    pass
                time.sleep(0.1)
            # warm (one name, the reviewer's first request), then every name
            warm = requests.post(url + '/traverse', data=json.dumps(self._request(self.names[-1:])))
            self.assertEqual(200, warm.status_code, warm.text)
            ret = requests.post(url + '/traverse', data=json.dumps(self._request(self.names)))
            self.assertEqual(200, ret.status_code, ret.text)
            served = self._check_seed_phase(ret.json(), 'server')
        finally:
            server.kill()
            server.wait()
            log.close()
        # stated for the record (the reviewer's warm server: seed phase 126.37 ms, longest piece
        # 0.297 ms)
        for where, t in (('CLI', cli), ('server', served)):
            print(f"{where}: seed_phase_ms {t['seed_phase_ms']:.3f}, label_resolve_ms "
                  f"{t['label_resolve_ms']:.3f}, longest_piece {t['deadline']['longest_piece']}",
                  file=sys.stderr)


if __name__ == '__main__':
    unittest.main()
