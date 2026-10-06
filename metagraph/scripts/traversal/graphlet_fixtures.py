#!/usr/bin/env python3
"""Build (or --check) the whole-document graphlet fixtures of DESIGN-traverse-graphlet.md
§2.7: api/python/tests/data/traverse/documents/<name>.{request.json,full.json,mgt,
graphlet.json}.

Two kinds of fixture live side by side:

  * CLI fixtures (CLI_FIXTURES below; request.json carries "fixture_index": {"cli": ...}).
    The script defines a handful of tiny indexes as FASTA built from seeded random blocks
    (INDEXES), builds each with the metagraph CLI (succinct graph, column_coord annotation,
    CoordToHeader sidecar), runs every request with detail full and with detail graphlet,
    and writes what the CLI printed: the .mgt is the graphlet body byte for byte, full.json
    and graphlet.json the two responses (timing off, so they are deterministic). Each
    fixture also names the record features it exists to show (`shows`); generating one that
    no longer shows them fails, so a walker change cannot silently turn a fixture into a
    different situation.
  * CLI comparison bodies (COMPARE_FIXTURES): body-only documents under documents/compare/
    named cli_*.mgt, pairs of retrievals of one seed whose comparisons (compare(), T39)
    the library tests check: a tuned and an exhaustive trie, and one locus at two radii
    and under a step cap.
  * Hand-made fixtures (HAND_MADE): situations the CLI cannot produce -- a Q record and an
    unknown counter (no resource stop exists yet), a resource-stopped delivery, and the
    clipped-merge case of §5.1 with its hand-derived claims. The .mgt body and full.json
    (written by hand from the walker's semantics) are the sources; the script derives their
    envelopes and graphlet.json.

usage:
    graphlet_fixtures.py [--dir DIR] [--check]                 # hand-made: derive / verify
    graphlet_fixtures.py --from-cli METAGRAPH [--check] [--only NAME ...] [--keep DIR]
                                                              # CLI: regenerate / verify
--check writes nothing and exits 1 when a file differs (a .mgt byte for byte, the JSON
files as parsed JSON). It does not import metagraph.traverse: the library is tested against
these files, not the other way round.
"""

import argparse
import copy
import json
import os
import random
import shutil
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DEFAULT_DIR = os.path.join(REPO, 'api', 'python', 'tests', 'data', 'traverse', 'documents')

WALK_RULE = ('a walk uses no (k+1)-mer edge twice on its own path, enters no seed node and '
             'takes no hairpin step (hand-made fixture: no index behind it)')

# the name every CLI fixture's index reports (H.index_ns, capabilities.index_ns). No manifest:
# its digest would be a digest of the serialised graph and annotation, which would make the
# fixtures depend on the serialisation format; index_fp is therefore null / '*'.
CLI_INDEX_NAME = 'fixtures'


# ===================================================================== CLI fixture indexes
#
# Every index is a few FASTA records assembled from random blocks (one seeded generator per
# index, so a block's bases never depend on another index). A record is a list of parts: a
# block name, a (block, start, end) slice, or a literal ('=', bases). Files are columns;
# record headers are the header labels (CoordToHeader).

INDEXES = {
    # the integration fixture's shape at k=15: acc1 = L1 E R1, acc2 = L1 E R2, acc3 = L2 E R1,
    # plus acc4, which carries only a prefix of E (a named label that does not support a seed
    # spanning E: dropped with its k-mer runs)
    'fork': {
        'k': 15, 'seed': 101,
        'blocks': {'L1': 40, 'L2': 40, 'E': 40, 'R1': 40, 'R2': 40, 'Z': 30},
        'files': {'refs.fa': [('acc1', ['L1', 'E', 'R1']),
                              ('acc2', ['L1', 'E', 'R2']),
                              ('acc3', ['L2', 'E', 'R1']),
                              ('acc4', [('E', 0, 25), 'Z'])]},
    },
    # two SNP bubbles to the right of the seed S, a third arm of the backbone to its left.
    # Alleles are chosen so that hB takes the higher child (the non-first parent) at both
    # bubbles: its run keeps the first merge's route stamp while its leaf entry carries the
    # second's (R.route_bp != T.route_bp). Column `both.fa` holds an hA and an hB copy, so
    # its lineage is ambiguous at each bubble and one of its two runs is closed by the merge
    # (an `m` end). The 25 bp spacers exceed k, so the bubbles are separate.
    'bubbles': {
        'k': 15, 'seed': 202,
        'blocks': {'W': 40, 'S': 30, 'P': 20, 'Q1': 25, 'Q2': 25, 'T': 40},
        'files': {'a.fa': [('hA', ['W', 'S', 'P', ('=', 'A'), 'Q1', ('=', 'C'), 'Q2', 'T'])],
                  'b.fa': [('hB', ['W', 'S', 'P', ('=', 'T'), 'Q1', ('=', 'G'), 'Q2', 'T'])],
                  'c.fa': [('hC', ['W', 'S', 'P', ('=', 'A'), 'Q1', ('=', 'G'), 'Q2', 'T'])],
                  'both.fa': [('hA2', ['W', 'S', 'P', ('=', 'A'), 'Q1', ('=', 'C'), 'Q2', 'T']),
                              ('hB2', ['W', 'S', 'P', ('=', 'T'), 'Q1', ('=', 'G'), 'Q2', 'T'])]},
    },
    # a switch chain A -> B -> {C, D}: A carries S T1, B the last 20 bp of T1 and T2, C and D
    # the last 20 bp of T2 and then T3 / T4. Walking right from S under A, A ends at the end
    # of T1 where only B goes on (a switch inside a segment); B ends at the end of T2, where
    # C and D go on along two successors (B's lineage on both: an ambiguous split whose
    # child sets {C} and {D} are disjoint, and whose switches sit on the children)
    'chain': {
        'k': 15, 'seed': 303,
        'blocks': {'S': 30, 'T1': 30, 'T2': 30, 'T3': 30, 'T4': 30},
        'files': {'chain.fa': [('A', ['S', 'T1']),
                               ('B', [('T1', 10, 30), 'T2']),
                               ('C', [('T2', 10, 30), 'T3']),
                               ('D', [('T2', 10, 30), 'T4'])]},
    },
    # A -> B -> A inside one segment: column pA holds S T1 and the last 20 bp of T2 with T3,
    # column pB the last 20 bp of T1 with T2. The walk leaves pA after T1, follows pB along
    # T2 and re-enters pA for T3, without a branch anywhere
    'reentry': {
        'k': 15, 'seed': 404,
        'blocks': {'S': 30, 'T1': 30, 'T2': 30, 'T3': 30},
        'files': {'pA.fa': [('a1', ['S', 'T1']), ('a2', [('T2', 10, 30), 'T3'])],
                  'pB.fa': [('b1', [('T1', 10, 30), 'T2'])]},
    },
    # two records named ACC1 in two files (two columns), both on the seed, so annotate mode
    # records them first and next to each other (the second L line front-codes the whole
    # name: an empty suffix); ACC2 leaves the first ACC1 20 bp into U1
    'same_name': {
        'k': 15, 'seed': 505,
        'blocks': {'S': 30, 'U1': 30, 'U2': 30, 'V': 20},
        'files': {'one.fa': [('ACC1', ['S', 'U1']), ('ACC2', [('U1', 0, 20), 'V'])],
                  'two.fa': [('ACC1', ['S', 'U2'])]},
    },
    # both stored-split counterexamples of §2.2 in one canonical locus. Right of S, a switch
    # chain over columns: a.fa -> b.fa at the end of T1, then b.fa ends at the end of T2
    # where c.fa and d.fa go on along two successors (b.fa's lineage on both: an AMBIGUOUS
    # split with DISJOINT child sets, its switches on the children). Left of S, a hairpin:
    # P = H8 + rc(H8) is a 16-mer equal to its own reverse complement, so at k = 15 the left
    # step that completes P turns onto the other strand. Both records of a.fa share rc(U) S
    # and the last 15 bases of P; A1 takes the palindromic step, A2 another base. a.fa is on
    # both predecessors, but a followed hairpin never counts toward ambiguity: a DIVERGENCE
    # with OVERLAPPING child sets
    'ambiguous_split': {
        'k': 15, 'seed': 707, 'mode': 'canonical',
        'blocks': {'S': 30, 'U': 10, 'H8': 8, 'V1': 30, 'V2': 30,
                   'T1': 30, 'T2': 30, 'T3': 30, 'T4': 30},
        'files': {'a.fa': [('A1', [('rc', 'V1'), 'H8', ('rc', 'H8'), ('rc', 'U'), 'S', 'T1']),
                           ('A2', [('rc', 'V2'), ('alt', 'H8', 0), ('H8', 1, 8), ('rc', 'H8'),
                                   ('rc', 'U'), 'S'])],
                  'b.fa': [('B', [('T1', 10, 30), 'T2'])],
                  'c.fa': [('C', [('T2', 10, 30), 'T3'])],
                  'd.fa': [('D', [('T2', 10, 30), 'T4'])]},
    },
    # quorum: four columns share S P; at the end of P, A carries x, y, z, B only x, C y and w.
    # min_successor_labels 2 refuses B (x's minority), so the split keeps A and C. x's lineage
    # was on two successors before the quorum, z's on one: when x and z both end at the end of
    # A (y goes on along D), their two label_lost ends look alike in every other field but
    # carry branches 1 and 0 -- the counterexample that made LabelRun::branches primary
    'quorum': {
        'k': 15, 'seed': 808,
        'blocks': {'S': 30, 'P': 10, 'A': 30, 'B': 30, 'C': 30, 'D': 30},
        'files': {'x.fa': [('x1', ['S', 'P', 'A']), ('x2', ['S', 'P', 'B'])],
                  'y.fa': [('y1', ['S', 'P', 'A', 'D']), ('y2', ['S', 'P', 'C'])],
                  'z.fa': [('z1', ['S', 'P', 'A'])],
                  'w.fa': [('w1', ['S', 'P', 'C'])]},
    },
    # caps and events in one locus: an indel bubble (hY carries an 8 bp insertion, so the
    # two routes reach the shared node at different depths: a revisit with a distance), a
    # repeat R that hX passes twice (the second entry would reuse R's edges: a blocked
    # successor, edge_reuse), and a loss-budget switch target hZ beyond the end of hX
    'caps': {
        'k': 15, 'seed': 606,
        'blocks': {'S': 30, 'P': 20, 'I': 8, 'Q': 30, 'R': 20, 'Y': 25, 'Z': 25, 'N': 30},
        'files': {'caps.fa': [('hX', ['S', 'P', 'Q', 'R', 'Y', 'R', 'Z']),
                              ('hY', ['S', 'P', 'I', 'Q', 'R', 'Z']),
                              ('hZ', [('Z', 5, 25), 'N'])]},
    },
    # a comb: a 260 bp backbone with a one-base dead-end tip at each of 230 positions (a
    # tip record is one k-mer: the backbone's last 14 bases there and another base), all
    # in one column. An exhaustive walk splits at every node, so the leaves' chains sum
    # quadratically: what detail full spells out and the memory budget charges (§14)
    'comb': {
        'k': 15, 'seed': 909,
        'blocks': {'B': 260},
        'files': {'comb.fa': [('backbone', ['B'])]
                  + [('tip%d' % i, [('B', i - 14, i), ('alt', 'B', i)]) for i in range(15, 245)]},
    },
    # a row-diff annotation with coordinates (stage 3 of DESIGN-traverse-graphlet.md §14.1,
    # the budget-aware reads): the path P (acc_path) and a record acc_rep that repeats the
    # k-mer P[51, 66) 120,000 times between random spacers. The k-mer before it on P carries
    # acc_path alone, but its row-diff successor carries 120,000 coordinates of acc_rep: a
    # row whose dependencies are dense while its result is tiny (§14 freeze gate)
    'dense': {
        'k': 15, 'seed': 1111, 'anno': 'row_diff_brwt_coord',
        'blocks': {'P': 120},
        'files': {'path.fa': [('acc_path', ['P'])],
                  'rep.fa': [('acc_rep', [('repeat', 'P', 51, 66, 120000, 10, 1112)])]},
    },
}


def blocks_of(spec):
    rng = random.Random(spec['seed'])
    # sorted, so adding a block appends to the generator's stream only for later names
    return {name: ''.join(rng.choice('ACGT') for _ in range(n))
            for name, n in sorted(spec['blocks'].items())}


def reverse_complement(s):
    return s[::-1].translate(str.maketrans('ACGT', 'TGCA'))


def assemble(parts, blocks):
    """A part is a block name, (block, start, end), ('=', literal bases), ('rc', block[,
    start, end]) -- a slice of the block's reverse complement -- or ('alt', block, i): one
    base other than the block's base i, so that a record leaves the others exactly there."""
    out = []
    for p in parts:
        if isinstance(p, str):
            out.append(blocks[p])
        elif p[0] == 'repeat':
            # ('repeat', block, start, end, times, spacer, seed): the slice |times| times,
            # each followed by |spacer| bases of a generator of its own
            _, block, start, end, times, spacer, seed = p
            rng = random.Random(seed)
            for _ in range(times):
                out.append(blocks[block][start:end])
                out.append(''.join(rng.choice('ACGT') for _ in range(spacer)))
        elif p[0] == '=':
            out.append(p[1])
        elif p[0] == 'rc':
            rc = reverse_complement(blocks[p[1]])
            out.append(rc[p[2]:p[3]] if len(p) > 2 else rc)
        elif p[0] == 'alt':
            out.append('ACGT'[('ACGT'.index(blocks[p[1]][p[2]]) + 1) % 4])
        else:
            out.append(blocks[p[0]][p[1]:p[2]])
    return ''.join(out)


def resolve_cli(path):
    """The CLI binary as an absolute path: build_index() runs it with cwd = the index
    directory, where a relative path ('build/metagraph') no longer resolves. A bare name
    is left to the PATH lookup."""
    if os.path.dirname(path):
        return os.path.abspath(path)
    return shutil.which(path) or path


def build_index(name, metagraph, root):
    """Build index |name| under root/name with the CLI; -> (graph, annotation) paths. The
    FASTA file names are passed relative (cwd = the index directory), so the column names,
    and with them index_meta_fp, do not depend on where the index is built."""
    spec = INDEXES[name]
    blocks = blocks_of(spec)
    d = os.path.join(root, name)
    os.makedirs(d, exist_ok=True)
    files = list(spec['files'])
    for fname, records in spec['files'].items():
        with open(os.path.join(d, fname), 'w') as f:
            for header, parts in records:
                f.write('>%s\n%s\n' % (header, assemble(parts, blocks)))
    k = str(spec['k'])
    if spec.get('anno') == 'row_diff_brwt_coord':
        # the row-diff pipeline (as integration_tests/base.py builds row_diff_brwt_coord):
        # the coordinate columns diffed along the graph in three stages, then compressed
        rd = [metagraph, 'transform_anno', '-p', '1', '--anno-type', 'row_diff', '--coordinates',
              '-o', 'annotation', '-i', 'graph.dbg', 'annotation.column.annodbg']
        transform = [rd, rd + ['--row-diff-stage', '1'], rd + ['--row-diff-stage', '2'],
                     [metagraph, 'transform_anno', '-p', '1', '--anno-type', 'row_diff_brwt_coord',
                      '--greedy', '-o', 'annotation', '-i', 'graph.dbg',
                      'annotation.column.annodbg'],
                     [metagraph, 'relax_brwt', '-p', '1', '-o', 'annotation',
                      'annotation.row_diff_brwt_coord.annodbg']]
    else:
        transform = [[metagraph, 'transform_anno', '-p', '1', '--anno-type', 'column_coord',
                      '--coordinates', '-o', 'annotation', 'annotation.column.annodbg']]
    for cmd in ([[metagraph, 'build', '-p', '1', '--mode', spec.get('mode', 'basic'),
                  '--graph', 'succinct', '-k', k, '-o', 'graph'] + files,
                 [metagraph, 'annotate', '-p', '1', '--anno-filename', '-i', 'graph.dbg',
                  '--anno-type', 'column', '--coordinates', '-o', 'annotation'] + files]
                + transform
                + [[metagraph, 'annotate', '-p', '1', '-i', 'graph.dbg', '--anno-filename',
                    '--index-header-coords', '-o', 'annotation'] + files]):
        p = subprocess.run(cmd, cwd=d, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if p.returncode:
            raise RuntimeError('%s failed:\n%s' % (' '.join(cmd), p.stderr.decode()[-2000:]))
    return os.path.join(d, 'graph.dbg'), annotation_of(name, root)


def annotation_of(index, root):
    """The annotation file of a built fixture index."""
    anno = INDEXES[index].get('anno', 'column_coord')
    return os.path.join(root, index, 'annotation.%s.annodbg' % anno)


# ===================================================================== CLI fixtures
#
# name -> (index, what it shows, the record features it must contain, request). A seed's
# sequence is given as parts (resolved against the index's blocks); timing is always off.
# Features (see features()): 'R:<end token>', 'E:<type>', 'G:merge', 'G:partition',
# 'G:split_amb' / 'G:split_div', 'G:no_bases', 'G:first_base', 'G:zero_length', 'T:<path
# reason>', 'T:extras', 'C', 'C:sequence', 'X', 'V:refused', 'P', 'P:cut', 'K:<kind>',
# 'L:empty_suffix', 'R:zero_length', 'R:route_split' (a run whose T.route_bp differs from its
# R.route_bp), 'E:revisit_distance', 'E:revisit_same', 'O:<4 codes>', 'A:no_bins' (an arm
# whose A record is followed by no B: it stopped before its first level).

CLI_FIXTURES = {
    'linear': {
        'index': 'fork',
        'doc': 'both arms, quorum 2: a minority end at each arm\'s fork with its V refusal '
               '(on the right at the seed boundary: a zero-length Rm run; L1 and L2 share their '
               'last base, so the left fork is at 1); acc4, named but carrying only part of '
               'the seed, is dropped (X)',
        'shows': ['R:Rm', 'R:zero_length', 'V:refused', 'X', 'R:D', 'O:ccci'],
        'request': {'seeds': [{'seed_id': 'linear', 'sequence': ['E'],
                               'labels': ['acc1', 'acc2', 'acc3', 'acc4']}],
                    'strategy': {'direction': 'both',
                                 'branching': {'min_successor_labels': 2},
                                 'bounds': {'max_extension_bp': 50},
                                 'output': {'profile_bin_bp': 20}}},
    },
    'fork': {
        'index': 'fork',
        'doc': 'both arms split (divergences); the radius (10) is below k (15), so both '
               'arms\' continuations cross into the seed',
        'shows': ['G:split_div', 'T:X', 'C', 'R:X'],
        'request': {'seeds': [{'seed_id': 'fork', 'sequence': ['E']}],
                    'strategy': {'direction': 'both', 'bounds': {'max_extension_bp': 10},
                                 'output': {'profile_bin_bp': 20}}},
    },
    'merge': {
        'index': 'bubbles',
        'doc': 'column labels through two SNP bubbles under on_reconverge merge: merge '
               'partitions (labels_via_parent), an ambiguous split (both.fa on both alleles) '
               'with an m closure at each merge; each merge\'s first parent is the one carried '
               'by the most labels (R21 (4)): at 36 the A allele (a.fa, c.fa, both.fa), at 62 '
               'the G allele (b.fa, c.fa, both.fa), which arrived second, so b.fa is routed in '
               'at 36 and a.fa at 62; continuations contained in the flank',
        'shows': ['G:merge', 'G:partition', 'G:split_amb', 'R:m', 'T:extras', 'C', 'O:pcci'],
        'request': {'seeds': [{'seed_id': 'merge', 'sequence': ['S']}],
                    'strategy': {'direction': 'right',
                                 'labels': {'seed_label_kind': 'column'},
                                 'branching': {'max_label_branches': 'unlimited'},
                                 'bounds': {'max_extension_bp': 100},
                                 'output': {'profile_bin_bp': 50, 'continuation_bp': 40}}},
    },
    'merge_ties': {
        'index': 'bubbles',
        'doc': 'the merge fixture without c.fa: both alleles of each bubble are carried by two '
               'labels, so each merge keeps its first-arrived parent first (a tie, R21 (4)), and '
               'b.fa goes through two non-first-parent merges (R.route_bp, the earliest stamp, '
               '!= T.route_bp, the latest)',
        'shows': ['G:merge', 'G:partition', 'G:split_amb', 'R:m', 'R:route_split', 'T:extras'],
        'request': {'seeds': [{'seed_id': 'merge_ties', 'sequence': ['S'],
                               'labels': ['a.fa', 'b.fa', 'both.fa']}],
                    'strategy': {'direction': 'right',
                                 'labels': {'seed_label_kind': 'column'},
                                 'branching': {'max_label_branches': 'unlimited'},
                                 'bounds': {'max_extension_bp': 100},
                                 'output': {'profile_bin_bp': 50, 'continuation_bp': 40}}},
    },
    'annotate': {
        'index': 'bubbles',
        'doc': 'annotate mode, merge: P runs with cut lists (cap 2), merges whose partition '
               'is one empty list per parent, contained continuations on both arms',
        'shows': ['P', 'P:cut', 'G:merge', 'K:label_lists', 'K:inexact_counts', 'T:X', 'C',
                  'O:pcli'],
        'request': {'seeds': [{'seed_id': 'annotate', 'sequence': ['S']}],
                    'strategy': {'direction': 'both',
                                 'labels': {'mode': 'annotate', 'max_labels_per_node': 2,
                                            'seed_label_kind': 'column'},
                                 'branching': {'on_reconverge': 'merge'},
                                 'bounds': {'max_extension_bp': 70},
                                 'output': {'profile_bin_bp': 50, 'continuation_bp': 20}}},
    },
    'switch_chain': {
        'index': 'chain',
        'doc': 'A -> B -> {C, D}: a silent Lw end inside a segment (switch on the same '
               'segment), and one at an ambiguous split with disjoint child sets whose '
               'switches sit on the children; both walks stop at the radius 20 bp after the '
               'last switch (continuations with loss_used 2)',
        'shows': ['R:Lw', 'E:switch', 'G:split_amb', 'R:X', 'T:extras', 'C'],
        'request': {'seeds': [{'seed_id': 'switch_chain', 'sequence': ['S'], 'labels': ['A']}],
                    'strategy': {'direction': 'right',
                                 'labels': {'extra': ['B', 'C', 'D'],
                                            'change_cost': {'model': 'constant', 'value': 1},
                                            'loss_budget': 3},
                                 'branching': {'max_label_branches': 1},
                                 'bounds': {'max_extension_bp': 80},
                                 'output': {'profile_bin_bp': 50}}},
    },
    'reentry': {
        'index': 'reentry',
        'doc': 'A -> B -> A inside one segment (the G.end counterexample to v1\'s set formula)',
        'shows': ['R:Lw', 'E:switch', 'R:D'],
        'request': {'seeds': [{'seed_id': 'reentry', 'sequence': ['S']}],
                    'strategy': {'direction': 'right',
                                 'labels': {'seed_label_kind': 'column', 'extra': ['pB.fa'],
                                            'change_cost': {'model': 'constant', 'value': 1},
                                            'loss_budget': 2},
                                 'bounds': {'max_extension_bp': 200},
                                 'output': {'profile_bin_bp': 50}}},
    },
    'same_name': {
        'index': 'same_name',
        'doc': 'annotate: two headers named ACC1 from different columns (an L line ending in '
               'a space)',
        'shows': ['L:empty_suffix', 'P', 'G:split_div'],
        'request': {'seeds': [{'seed_id': 'same_name', 'sequence': ['S']}],
                    'strategy': {'direction': 'right',
                                 'labels': {'mode': 'annotate', 'seed_label_kind': 'header'},
                                 'bounds': {'max_extension_bp': 40},
                                 'output': {'profile_bin_bp': 50}}},
    },
    'ambiguous_split': {
        'index': 'ambiguous_split',
        'doc': 'canonical graph, column labels: an ambiguous split with disjoint child sets '
               'after a switch (right); a followed hairpin at a divergence whose child sets '
               'overlap (left)',
        'shows': ['G:split_amb', 'G:split_div', 'E:hairpin_followed', 'E:switch', 'R:Lw'],
        'request': {'seeds': [{'seed_id': 'ambiguous_split', 'sequence': ['S']}],
                    'strategy': {'direction': 'both',
                                 'labels': {'seed_label_kind': 'column',
                                            'extra': ['b.fa', 'c.fa', 'd.fa'],
                                            'change_cost': {'model': 'constant', 'value': 1},
                                            'loss_budget': 2},
                                 'branching': {'hairpins': 'follow', 'max_label_branches': 1},
                                 'bounds': {'max_extension_bp': 100},
                                 'output': {'profile_bin_bp': 50}}},
    },
    'quorum': {
        'index': 'quorum',
        'doc': 'quorum 2 at a three-way fork: the minority successor is refused (V) and the '
               'split keeps two children; two interior label_lost ends at one position whose '
               'branch counts (1 and 0) differ only through the refused successor',
        'shows': ['V:refused', 'G:split_amb', 'R:L', 'R:branch_counts_differ'],
        'request': {'seeds': [{'seed_id': 'quorum', 'sequence': ['S']}],
                    'strategy': {'direction': 'right',
                                 'labels': {'seed_label_kind': 'column'},
                                 'branching': {'min_successor_labels': 2,
                                               'max_label_branches': 1},
                                 'bounds': {'max_extension_bp': 100},
                                 'output': {'profile_bin_bp': 50}}},
    },
    'caps': {
        'index': 'caps',
        'doc': 'caps and events: a max_steps truncation, branch events cut at 1, sequences '
               'false (first bases, explicit continuation sequences), a loss_budget end, a '
               'blocked edge_reuse successor and a revisit at another distance',
        'shows': ['K:walk_domain', 'K:branch_events', 'G:no_bases', 'C:sequence', 'R:B',
                  'E:blocked', 'E:revisit_distance', 'T:S', 'O:pxci'],
        'request': {'seeds': [{'seed_id': 'caps', 'sequence': ['S'], 'labels': ['hX', 'hY']}],
                    'strategy': {'direction': 'right',
                                 'labels': {'extra': ['hZ'],
                                            'change_cost': {'model': 'table',
                                                            'entries': [['hX', 'hZ', 2],
                                                                        ['hY', 'hZ', 0.5]],
                                                            'default': 'forbid'},
                                            'loss_budget': 1},
                                 'branching': {'on_reconverge': 'keep',
                                               'max_label_branches': 'unlimited'},
                                 'bounds': {'max_extension_bp': 200, 'max_steps': 240},
                                 'output': {'sequences': False, 'max_branch_events': 1,
                                            'profile_bin_bp': 50}}},
    },
    'no_bins': {
        'index': 'bubbles',
        'doc': 'annotate, both arms, max_steps 1 below the left boundary\'s two successors (the '
               'seed Q1 sits between two bubbles): the cap trips on the left arm at depth 0 and '
               'censors the right root before it ran a level, so the right arm has no growth '
               'bin (A B* V*: zero B records are valid)',
        'shows': ['A:no_bins', 'T:S', 'G:zero_length', 'K:walk_domain', 'O:pcci'],
        'request': {'seeds': [{'seed_id': 'no_bins', 'sequence': ['Q1']}],
                    'strategy': {'direction': 'both',
                                 'labels': {'mode': 'annotate', 'seed_label_kind': 'column'},
                                 'bounds': {'max_extension_bp': 50, 'max_steps': 1},
                                 'output': {'profile_bin_bp': 20}}},
    },
}

# body-only CLI retrievals for the comparison tests: documents/compare/<name>.mgt. The
# quorum pair is the T39 prefix-subset case (the tuned trie refuses the minority at two
# forks and ends two labels inside a walk that goes on with others); the bubbles triple
# compares one locus at radius 60 and 100 and under a step cap (claims ending at the
# comparison depth). on_reconverge keep everywhere, so the scopes agree.
_QUORUM = {'seeds': [{'seed_id': 'q', 'sequence': ['S']}],
           'strategy': {'direction': 'right', 'labels': {'seed_label_kind': 'column'},
                        'bounds': {'max_extension_bp': 100}}}
_BUBBLES = {'seeds': [{'seed_id': 'b', 'sequence': ['S']}],
            'strategy': {'direction': 'right', 'labels': {'seed_label_kind': 'column'},
                         'branching': {'max_label_branches': 'unlimited',
                                       'on_reconverge': 'keep'},
                         'bounds': {'max_extension_bp': 100}}}


def _variant(base, **blocks):
    r = copy.deepcopy(base)
    for block, values in blocks.items():
        if isinstance(values, dict):
            r['strategy'].setdefault(block, {}).update(values)
        else:
            r['strategy'][block] = values
    return r


COMPARE_FIXTURES = {
    'cli_quorum_tuned': ('quorum', _variant(_QUORUM, branching={
        'min_successor_labels': 2, 'max_label_branches': 1, 'on_reconverge': 'keep'})),
    'cli_quorum_exhaustive': ('quorum', _variant(_QUORUM, exhaustive=True)),
    'cli_bubbles_r60': ('bubbles', _variant(_BUBBLES, bounds={'max_extension_bp': 60})),
    'cli_bubbles_r100': ('bubbles', _BUBBLES),
    'cli_bubbles_steps150': ('bubbles', _variant(_BUBBLES, bounds={'max_steps': 150})),
}

# hand-made: what the CLI cannot produce today
HAND_MADE = ('clipped_merge', 'limits')

# CLI fixtures of the request budgets (DESIGN-traverse-graphlet.md §14, stage 2), written to
# documents/budgets/ (the top-level fixture list stays as it is). A work stop is the same
# walk in every detail, so its full and graphlet responses are both written; a memory
# stop is charged in the requested detail (the output is part of the budget), so the same
# request stops at another depth in another detail and each memory stop is written in its
# own detail only.
BUDGET_DIR = 'budgets'
BUDGET_FIXTURES = {
    'work_stop': {
        'index': 'bubbles',
        'doc': 'the merge fixture under bounds.max_work_units 1500: the walk stops between two '
               'heads at the work budget (Q work, resource_limit leaves and runs, a walk_domain '
               'limitation on the budget, the work_units counter)',
        'shows': ['Q:work', 'T:Y', 'R:Y', 'K:walk_domain'],
        'details': ('full', 'graphlet'),
        'request': {'seeds': [{'seed_id': 'work_stop', 'sequence': ['S']}],
                    'strategy': {'direction': 'right',
                                 'labels': {'seed_label_kind': 'column'},
                                 'branching': {'max_label_branches': 'unlimited'},
                                 'bounds': {'max_extension_bp': 100, 'max_work_units': 1500},
                                 'output': {'profile_bin_bp': 50, 'continuation_bp': 40}}},
    },
    'memory_soft': {
        'index': 'quorum',
        'doc': 'the quorum fixture under bounds.max_memory_mb 1, which admits every head: no '
               'stop, but the seed-level memory_bound_soft limitation (annotation decoding is '
               'not charged yet)',
        'shows': ['K:memory_bound_soft', 'V:refused'],
        'details': ('full', 'graphlet'),
        'request': {'seeds': [{'seed_id': 'memory_soft', 'sequence': ['S']}],
                    'strategy': {'direction': 'right',
                                 'labels': {'seed_label_kind': 'column'},
                                 'branching': {'min_successor_labels': 2,
                                               'max_label_branches': 1},
                                 'bounds': {'max_extension_bp': 100, 'max_memory_mb': 1},
                                 'output': {'profile_bin_bp': 50}}},
    },
    'memory_stop': {
        'index': 'comb',
        'doc': 'the exhaustive comb under bounds.max_memory_mb 1, detail graphlet: the budget '
               'stops the walk (Q memory, walk_domain on the budget, memory_bound_soft)',
        'shows': ['Q:memory', 'T:Y', 'K:walk_domain', 'K:memory_bound_soft'],
        'details': ('graphlet',),
        'request': {'seeds': [{'seed_id': 'memory_stop', 'sequence': [('B', 0, 15)]}],
                    'strategy': {'direction': 'right', 'exhaustive': True,
                                 'labels': {'seed_label_kind': 'column'},
                                 'bounds': {'max_extension_bp': 240, 'max_memory_mb': 1},
                                 'output': {'profile_bin_bp': 50, 'continuation_bp': 0,
                                            'max_branch_events': 'unlimited'}}},
    },
    'decode_stop': {
        'index': 'dense',
        'doc': 'the dense row-diff index (budget-aware annotation reads, stage 3) under '
               'bounds.max_memory_mb 1: the level that reads the row before the repeated '
               'k-mer cannot hold its row-diff dependency rows, so the walk stops there in '
               'phase annotation_decode (Q memory annotation_decode, its levers, a walk_domain '
               'on the budget, memory_bound_soft naming only what is still uncharged)',
        'shows': ['Q:memory', 'Q:annotation_decode', 'T:Y', 'K:walk_domain',
                  'K:memory_bound_soft'],
        'details': ('graphlet',),
        'request': {'seeds': [{'seed_id': 'decode_stop', 'sequence': [('P', 0, 30)],
                               'labels': ['acc_path']}],
                    'strategy': {'direction': 'right',
                                 'bounds': {'max_extension_bp': 60, 'max_memory_mb': 1},
                                 'output': {'profile_bin_bp': 50}}},
    },
    'memory_stop_full': {
        'index': 'comb',
        'doc': 'the memory_stop request in detail full: the JSON spells every leaf\'s chain '
               '(quadratic on a comb), so the same budget stops it earlier than the graphlet',
        'shows': [],
        'details': ('full',),
        'request': {'seeds': [{'seed_id': 'memory_stop', 'sequence': [('B', 0, 15)]}],
                    'strategy': {'direction': 'right', 'exhaustive': True,
                                 'labels': {'seed_label_kind': 'column'},
                                 'bounds': {'max_extension_bp': 240, 'max_memory_mb': 1},
                                 'output': {'profile_bin_bp': 50, 'continuation_bp': 0,
                                            'max_branch_events': 'unlimited'}}},
    },
}


def materialize(request, blocks):
    """The request as sent: seed sequences resolved, timing off (deterministic JSON)."""
    r = copy.deepcopy(request)
    for s in r['seeds']:
        s['sequence'] = assemble(s['sequence'], blocks)
    r['strategy'].setdefault('output', {})['timing'] = False
    return r


def features(body):
    """The record features of one MGT body (a cheap reading of the text, no library)."""
    f = set()
    lines = body.split('\n')
    t_route = {}     # (arm, leaf segment) -> {label: route_bp} from T extras
    runs = []        # (arm, segment, label, route_bp)
    ends = {}        # (arm, segment, to_bp, end token, from_bp) -> {branches} of ended runs
    arm = None
    seg = -1
    run = -1
    after_a = False  # the previous record was an A: is the arm's first record a B?
    for line in lines:
        if not line:
            continue
        x = line.split(' ')
        tag = x[0]
        if after_a and tag != 'B':
            f.add('A:no_bins')
        after_a = tag == 'A'
        if tag == 'A':
            arm, seg, run = x[1], -1, -1
        elif tag == 'O':
            f.add('O:' + ''.join(x[1:5]))
        elif tag == 'Q':
            f.add('Q:' + x[2])
            f.add('Q:' + x[3])
        elif tag == 'K':
            f.add('K:' + x[2])
        elif tag == 'X':
            f.add('X')
        elif tag == 'L' and line.endswith(' ') and x[4] != '0':
            f.add('L:empty_suffix')
        elif tag == 'V' and x[7] != '.':
            f.add('V:refused')
        elif tag == 'G':
            seg += 1
            if x[1] != '*' and ',' in x[1]:
                f.add('G:merge')
            if x[7] != '*':
                f.add('G:partition')
            if x[8] == '1':
                f.add('G:split_amb')
            elif x[8] == '0':
                f.add('G:split_div')
            if x[9] != '*':
                f.add('G:first_base')
            if x[10] == '*':
                f.add('G:no_bases')
            if x[3] == '0':
                f.add('G:zero_length')
        elif tag == 'P':
            f.add('P')
            if int(x[3]) > len(x[4].split(',')) and x[4] != '!':
                f.add('P:cut')
        elif tag == 'E':
            f.add('E:' + {'s': 'switch', 'b': 'blocked', 'h': 'hairpin', 'v': 'revisit',
                          't': 'tip', 'u': 'bubble'}[x[2]])
            if x[2] == 'v':
                f.add('E:revisit_same' if x[4] == '=' else 'E:revisit_distance')
            if x[2] == 'h':
                f.add('E:hairpin_followed' if x[-1] == 'f' else 'E:hairpin_skipped')
        elif tag == 'T':
            f.add('T:' + x[1])
            if x[2] != '.':
                f.add('T:extras')
                for item in x[2].split(','):
                    lab, _, _, route = item.split(':')
                    t_route[(arm, seg, int(lab))] = int(route)
        elif tag == 'C':
            f.add('C:sequence' if len(x) > 5 else 'C')
        elif tag == 'R':
            run += 1
            f.add('R:' + x[5])
            if x[3] == x[4]:
                f.add('R:zero_length')
            runs.append((arm, int(x[1]), int(x[2]), int(x[6])))
            if x[5] not in ('Lw', 'm'):
                ends.setdefault((arm, x[1], x[4], x[5], x[3]), set()).add(x[10])
    for a, s, label, route in runs:
        if (a, s, label) in t_route and t_route[(a, s, label)] != route:
            f.add('R:route_split')
    for a, x in ends.items():
        if len(x) > 1:
            f.add('R:branch_counts_differ')
    return f


# ===================================================================== hand-made envelopes

def normalized_strategy(strategy, detail):
    """strategy_to_json() of src/cli/traverse.cpp over the request's strategy: every
    default filled in, the mode/preset dependent ones as parse_traverse_request() sets
    them, plus output.detail/timing (v2 echo) and clamped."""
    st = copy.deepcopy(strategy or {})
    labels = st.get('labels', {})
    mode = labels.get('mode', 'constrain')
    exhaustive = st.get('exhaustive', False)
    annotate = mode == 'annotate'
    max_seed = labels.get('max_seed_labels', 1000)
    out = {
        'schema_version': 1, 'direction': st.get('direction', 'both'),
        'support': st.get('support', 'kmer'), 'exhaustive': exhaustive,
        'labels': {
            'mode': mode,
            'max_labels_per_node': labels.get('max_labels_per_node',
                                              max(64, max_seed) if annotate else 64),
            'extra': labels.get('extra', []),
            'change_cost': labels.get('change_cost', {'model': 'forbid'}),
            'loss_budget': float(labels.get('loss_budget', 0)),
            'switch_on': labels.get('switch_on', 'loss'),
            'max_switch_sources': labels.get('max_switch_sources',
                                             'unlimited' if exhaustive else 64),
            'max_seed_labels': max_seed,
            'seed_label_kind': labels.get('seed_label_kind', 'auto'),
        },
        'branching': {
            'max_label_branches': 'unlimited' if exhaustive or annotate else 0,
            'min_successor_labels': 1, 'min_successor_fraction': 0.0,
            'tip_window_bp': 0, 'bubble_window_bp': 0,
            'on_reconverge': 'keep' if exhaustive or annotate else 'merge',
            'hairpins': 'skip',
            'max_splits_per_path': 'unlimited' if exhaustive or annotate else 64,
        },
        'bounds': {'max_extension_bp': 5000, 'min_live_labels': 1, 'max_steps': 1000000,
                   'max_live_paths': 1000, 'max_paths': 10000, 'max_output_bp': 2000000,
                   'time_budget_ms': 30000.0},
        'frontier': {'order': 'breadth_first', 'on_overflow': 'stop'},
        'output': {'detail': detail, 'sequences': True, 'profile_bin_bp': 100,
                   'max_branch_events': 100, 'continuation_bp': 1000, 'timing': True},
        'annotation': {'batch_kmers': 64},
    }
    for block in ('branching', 'bounds', 'frontier', 'output', 'annotation'):
        for k, v in st.get(block, {}).items():
            out[block][k] = v
    out['output']['detail'] = detail
    out['bounds']['time_budget_ms'] = float(out['bounds']['time_budget_ms'])
    out['branching']['min_successor_fraction'] = float(
        out['branching']['min_successor_fraction'])
    out['clamped'] = []
    return out


def envelope(request, detail, index):
    caps = {
        'schema_version': 1, 'k': index['k'], 'regime': index.get('regime', 'basic'),
        'alphabet': index.get('alphabet', '$ACGT'), 'num_labels': index.get('num_labels', 1),
        'has_coordinates': True, 'has_coord_to_header': True, 'supports_trace': True,
        'cost_models_available': ['forbid', 'constant', 'table'],
        'label_modes': ['constrain', 'annotate'], 'direct_access': False,
        'release': index.get('release', ''),
        'graphlet_format': 1, 'detail_levels': ['summary', 'tree', 'full', 'graphlet'],
        'index_ns': index.get('ns'), 'index_fp': index.get('fp'),
        'index_meta_fp': index.get('meta_fp'),
    }
    return {'release': index.get('release', ''), 'capabilities': caps,
            'strategy': normalized_strategy(request.get('strategy'), detail),
            'walk_rule': WALK_RULE, 'algorithm_version': 'traverse-0.2'}


def arm_summary(arm):
    segs = arm['segments']
    leaves_by_reason = {}
    for p in arm['paths']:
        r = p.get('path_reason', 'semantic')
        leaves_by_reason[r] = leaves_by_reason.get(r, 0) + 1
    label_ends = {}
    for b in arm['growth']:
        for r, n in (b.get('label_ends') or {}).items():
            label_ends[r] = label_ends.get(r, 0) + n
    j = {k: copy.deepcopy(arm[k]) for k in ('status', 'complete_to_bp', 'completeness_scope',
                                            'frontier_remaining', 'labels_per_node',
                                            'counters', 'evidence', 'limitations')}
    if 'cap_trigger' in arm:
        j['cap_trigger'] = copy.deepcopy(arm['cap_trigger'])
    j['branch_events_total'] = len(arm['branch_events']) + arm['branch_events_truncated']
    j['counts'] = {
        'segments': len(segs), 'leaves': len(arm['paths']), 'splits': len(arm['splits']),
        'merges': sum(1 for s in segs if len(s['parents']) > 1), 'runs': len(arm['runs']),
        'bases': sum(s['length_bp'] for s in segs),
        'max_bp': max((p['length_bp'] for p in arm['paths']), default=0),
        'label_ends': label_ends, 'leaves_by_reason': leaves_by_reason}
    return j


def seed_summary(result, body):
    """§3: no names, no per-leaf rows, no label_summary, no growth, no label_dict."""
    seed = result['seed']
    out = {'seed': {k: seed[k] for k in ('seed_id', 'validated_seed_id', 'seed_id_mismatch',
                                          'length_bp', 'num_kmers', 'labels_from_seed',
                                          'labels_supporting_total', 'labels_dropped',
                                          'labels_dropped_digest')},
           'label_mode': result['label_mode']}
    out['seed']['num_labels'] = len(result['label_dict'])
    out['seed']['num_seed_labels'] = len(seed['labels'])
    if result.get('duplicate'):
        out['duplicate'] = True
    out['arms'] = {side: arm_summary(arm) for side, arm in result['arms'].items()}
    out['annotation'] = copy.deepcopy(result['annotation'])
    if 'timing' in result:
        out['timing'] = copy.deepcopy(result['timing'])
    out['outcome'] = copy.deepcopy(result['outcome'])
    if 'resource_stop' in result:
        out['resource_stop'] = copy.deepcopy(result['resource_stop'])
    out['limitations'] = copy.deepcopy(result['limitations'])
    out['graphlet'] = body
    out['graphlet_bytes'] = len(body.encode('utf-8'))
    out['graphlet_lines'] = body.count('\n')
    return out


# ===================================================================== files

def read_json(path):
    with open(path, 'r', encoding='utf-8') as f:
        return json.load(f)


def dumps(obj):
    return json.dumps(obj, indent=1, ensure_ascii=False, sort_keys=False) + '\n'


def read_text(path):
    with open(path, 'r', encoding='utf-8', newline='') as f:
        return f.read()


def build_hand_made(d, name):
    """-> {filename: text} that this hand-made fixture's derived files must hold."""
    request = read_json(os.path.join(d, name + '.request.json'))
    full = read_json(os.path.join(d, name + '.full.json'))
    body = read_text(os.path.join(d, name + '.mgt'))
    index = request.get('fixture_index') or {}
    graphlet = dict(envelope(request, 'graphlet', index))
    graphlet['results'] = [seed_summary(full['results'][0], body)]
    full_with_env = dict(envelope(request, 'full', index))
    full_with_env['results'] = full['results']
    return {name + '.graphlet.json': dumps(graphlet), name + '.full.json': dumps(full_with_env)}


def run_traverse(metagraph, graph, annotation, request):
    with tempfile.NamedTemporaryFile('w', suffix='.json', delete=False) as f:
        json.dump(request, f)
    try:
        p = subprocess.run([metagraph, 'traverse', '--index-name', CLI_INDEX_NAME,
                            '-i', graph, '-a', annotation, f.name],
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    finally:
        os.unlink(f.name)
    if p.returncode:
        raise RuntimeError('traverse failed (%d): %s' % (p.returncode,
                                                          p.stdout.decode()[-1000:]
                                                          + p.stderr.decode()[-1000:]))
    return json.loads(p.stdout.decode('utf-8'))


def build_cli(metagraph, root, name):
    """-> {filename: text} for CLI fixture |name|, and the features its body shows."""
    fx = CLI_FIXTURES[name]
    spec = INDEXES[fx['index']]
    graph = os.path.join(root, fx['index'], 'graph.dbg')
    if not os.path.exists(graph):
        build_index(fx['index'], metagraph, root)
    annotation = annotation_of(fx['index'], root)
    request = materialize(fx['request'], blocks_of(spec))
    outs = {}
    for detail in ('full', 'graphlet'):
        r = copy.deepcopy(request)
        r['strategy']['output']['detail'] = detail
        outs[detail] = run_traverse(metagraph, graph, annotation, r)
    body = outs['graphlet']['results'][0]['graphlet']
    stored = dict(request, fixture_index={'cli': fx['index'], 'k': spec['k'], 'doc': fx['doc']})
    files = {name + '.request.json': dumps(stored), name + '.mgt': body,
             name + '.full.json': dumps(outs['full']),
             name + '.graphlet.json': dumps(outs['graphlet'])}
    return files, features(body)


def build_budget(metagraph, root, name):
    """-> {filename: text} for budget fixture |name| (in BUDGET_DIR), and the features its
    body shows (none without a graphlet)."""
    fx = BUDGET_FIXTURES[name]
    spec = INDEXES[fx['index']]
    graph = os.path.join(root, fx['index'], 'graph.dbg')
    if not os.path.exists(graph):
        build_index(fx['index'], metagraph, root)
    annotation = annotation_of(fx['index'], root)
    request = materialize(fx['request'], blocks_of(spec))
    stored = dict(request, fixture_index={'cli': fx['index'], 'k': spec['k'], 'doc': fx['doc']})
    files = {name + '.request.json': dumps(stored)}
    shown = set()
    for detail in fx['details']:
        r = copy.deepcopy(request)
        r['strategy']['output']['detail'] = detail
        out = run_traverse(metagraph, graph, annotation, r)
        files['%s.%s.json' % (name, detail)] = dumps(out)
        if detail == 'graphlet':
            body = out['results'][0]['graphlet']
            files[name + '.mgt'] = body
            shown = features(body)
    return {os.path.join(BUDGET_DIR, f): text for f, text in files.items()}, shown


def build_compare(metagraph, root, name):
    """-> the body of comparison fixture |name| (detail graphlet, timing off)."""
    index, request = COMPARE_FIXTURES[name]
    graph = os.path.join(root, index, 'graph.dbg')
    if not os.path.exists(graph):
        build_index(index, metagraph, root)
    annotation = annotation_of(index, root)
    r = materialize(request, blocks_of(INDEXES[index]))
    r['strategy']['output']['detail'] = 'graphlet'
    return run_traverse(metagraph, graph, annotation, r)['results'][0]['graphlet']


def differs(path, text):
    if not os.path.exists(path):
        return True
    have = read_text(path)
    if have == text:
        return False
    if path.endswith('.json'):
        return json.loads(have) != json.loads(text)
    return True


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--dir', default=DEFAULT_DIR)
    ap.add_argument('--check', action='store_true', help='verify, write nothing')
    ap.add_argument('--from-cli', metavar='METAGRAPH',
                    help='regenerate (or with --check verify) the CLI fixtures')
    ap.add_argument('--only', nargs='*')
    ap.add_argument('--keep', metavar='DIR', help='build the indexes here and keep them')
    args = ap.parse_args(argv)
    stale = []
    if args.from_cli:
        args.from_cli = resolve_cli(args.from_cli)
        root = args.keep or tempfile.mkdtemp(prefix='graphlet-fixtures-')
        try:
            for name in CLI_FIXTURES:
                if args.only and name not in args.only:
                    continue
                files, shown = build_cli(args.from_cli, root, name)
                missing = sorted(set(CLI_FIXTURES[name]['shows']) - shown)
                if missing:
                    stale.append('%s no longer shows %s (shows %s)'
                                 % (name, missing, sorted(shown)))
                    continue
                for fname, text in files.items():
                    path = os.path.join(args.dir, fname)
                    if not differs(path, text):
                        continue
                    if args.check:
                        stale.append(fname)
                    else:
                        with open(path, 'w', encoding='utf-8', newline='') as f:
                            f.write(text)
                        print('wrote', fname)
            for name in BUDGET_FIXTURES:
                if args.only and name not in args.only:
                    continue
                files, shown = build_budget(args.from_cli, root, name)
                missing = sorted(set(BUDGET_FIXTURES[name]['shows']) - shown)
                if missing:
                    stale.append('%s no longer shows %s (shows %s)'
                                 % (name, missing, sorted(shown)))
                    continue
                os.makedirs(os.path.join(args.dir, BUDGET_DIR), exist_ok=True)
                for fname, text in files.items():
                    path = os.path.join(args.dir, fname)
                    if not differs(path, text):
                        continue
                    if args.check:
                        stale.append(fname)
                    else:
                        with open(path, 'w', encoding='utf-8', newline='') as f:
                            f.write(text)
                        print('wrote', fname)
            for name in COMPARE_FIXTURES:
                if args.only and name not in args.only:
                    continue
                fname = os.path.join('compare', name + '.mgt')
                text = build_compare(args.from_cli, root, name)
                path = os.path.join(args.dir, fname)
                if not differs(path, text):
                    continue
                if args.check:
                    stale.append(fname)
                else:
                    with open(path, 'w', encoding='utf-8', newline='') as f:
                        f.write(text)
                    print('wrote', fname)
        finally:
            if not args.keep:
                shutil.rmtree(root, ignore_errors=True)
    else:
        for name in HAND_MADE:
            if args.only and name not in args.only:
                continue
            for fname, text in build_hand_made(args.dir, name).items():
                path = os.path.join(args.dir, fname)
                if not differs(path, text):
                    continue
                if args.check:
                    stale.append(fname)
                else:
                    with open(path, 'w', encoding='utf-8', newline='') as f:
                        f.write(text)
                    print('wrote', fname)
    if stale:
        print('out of date:\n  ' + '\n  '.join(stale), file=sys.stderr)
        return 1
    if args.check:
        print('fixtures up to date (%s)' % ('CLI' if args.from_cli else 'hand-made'))
    return 0


if __name__ == '__main__':
    sys.exit(main())
