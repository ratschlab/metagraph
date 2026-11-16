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
import difflib
import json
import os
import platform
import random
import re
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
    # carry branches 1 and 0 -- the case that needs LabelRun::branches in the primary order
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
    # a row-diff annotation with coordinates (the budget-aware reads of
    # DESIGN-traverse-graphlet.md §14.1): the path P (acc_path) and a record acc_rep that
    # repeats the k-mer P[51, 66) 120,000 times between random spacers. The k-mer before it on
    # P carries acc_path alone, but its row-diff successor carries 120,000 coordinates of
    # acc_rep: a row whose dependencies are dense while its result is tiny (§14 freeze gate)
    'dense': {
        'k': 15, 'seed': 1111, 'anno': 'row_diff_brwt_coord',
        'blocks': {'P': 120},
        'files': {'path.fa': [('acc_path', ['P'])],
                  'rep.fa': [('acc_rep', [('repeat', 'P', 51, 66, 120000, 10, 1112)])]},
    },
    # ---- the stitching triples' indexes (STITCH_FIXTURES; DESIGN-graphlet-merge.md E.2 T4).
    # Each seed was chosen so that the records hold exactly the repeats named (no other
    # k-mer or (k+1)-mer repeats).
    # the external review's counterexample (Appendix G.4): one record S P U C U Y whose 15-mer
    # U is the one repeat, so the node u has two successors, C[0] (the loop) and Y[0]
    'loop': {
        'k': 15, 'seed': 1301,
        'blocks': {'S': 30, 'P': 40, 'U': 15, 'C': 45, 'Y': 30},
        'files': {'a.fa': [('r1', ['S', 'P', 'U', 'C', 'U', 'Y'])]},
    },
    # a repeat of one (k+1)-mer R preceded by two different bases: S P R Q R T. A walk along Q
    # enters R's first k-mer again by a new edge and then needs R's own edge (edge_reuse); a
    # child seeded with the walk blocks the new edge as rejoined_seed one base earlier
    'repeat': {
        'k': 15, 'seed': 1405,
        'blocks': {'S': 30, 'P': 30, 'R': 16, 'Q': 40, 'T': 20},
        'files': {'a.fa': [('r1', ['S', 'P', 'R', 'Q', 'R', 'T'])]},
    },
    # one label on three records: S M A1, S M A2 and S M A1[0:20] B, so the label is ambiguous
    # at the end of M and again 20 bp into A1 (branch events at two depths)
    'branchy': {
        'k': 15, 'seed': 1502,
        'blocks': {'S': 30, 'M': 40, 'A1': 60, 'A2': 10, 'B': 30},
        'files': {'a.fa': [('a1', ['S', 'M', 'A1']), ('a2', ['S', 'M', 'A2']),
                           ('a3', ['S', 'M', ('A1', 0, 20), 'B'])]},
    },
    # a.fa on both successors of the seed (S X1, S X2: an ambiguous split, X2 a 10 bp dead
    # end), b.fa on the X1 side only and shorter (45 bp)
    'whitelist': {
        'k': 15, 'seed': 1601,
        'blocks': {'S': 30, 'X1': 60, 'X2': 10},
        'files': {'a.fa': [('a1', ['S', 'X1']), ('a2', ['S', 'X2'])],
                  'b.fa': [('b1', ['S', ('X1', 0, 45)])]},
    },
    # x.fa on S T, z.fa on S T[0:20], y.fa on T[5:60]: z.fa is lost at 20 and switches to
    # y.fa, which is alive 5 bp later without covering a continuation's 15 bp tail
    'switch_tail': {
        'k': 15, 'seed': 1701,
        'blocks': {'S': 30, 'T': 60},
        'files': {'x.fa': [('x1', ['S', 'T'])], 'z.fa': [('z1', ['S', ('T', 0, 20)])],
                  'y.fa': [('y1', [('T', 5, 60)])]},
    },
    # the full Fable review's split-count counterexample (Appendix G.5, `splits45`), k 11:
    # A = S U VA R TA, B = S U VB, C = D = S WC R TC: the walk of A and C reconverges on R
    # through a merge whose parents hold one and two splits
    'splits': {
        'k': 11, 'seed': 1801,
        'blocks': {'S': 30, 'U': 10, 'VA': 10, 'VB': 8, 'WC': 20, 'R': 30, 'TA': 40, 'TC': 40},
        'files': {'a.fa': [('A', ['S', 'U', 'VA', 'R', 'TA'])],
                  'b.fa': [('B', ['S', 'U', 'VB'])],
                  'c.fa': [('C', ['S', 'WC', 'R', 'TC'])],
                  'd.fa': [('D', ['S', 'WC', 'R', 'TC'])]},
    },
    # the GPT review's circular record (Appendix G.6): its last 15 bases are its first, so
    # the continuation of a walk from either 15 bp seed is the other's seed
    'cycle': {
        'k': 15, 'seed': 1901, 'blocks': {},
        'files': {'A.fa': [('record', [('=', 'GGCCCCGCCGCAACTGACATCAAGTTGATTGGTCGTTAATGCGCACC'
                                              'TCTGGCCCCGCCGCAACT')])]},
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

# The stitching triples (DESIGN-graphlet-merge.md §3, Appendix E.2 T4), written to
# documents/stitch/: a parent at radius R, its child -- the continuation of one walk, run
# with the request parent.next_request(arm, [walk], bp=E, ...) builds (stitch_child_request()
# below follows its rules; test_traverse_stitch pins the two equal) -- and the uninterrupted
# retrieval at R + E. Some add a sibling (the continuation of another walk of the parent)
# or a grandchild (a continuation of the child's walk). Files: <name>.request.json (the
# requests, the child's among them) and <name>.<role>.json, the detail graphlet response of
# each role. Each fixture states what its walks must look like (expect: the continued
# walk's length L and continuation length n, the walk lengths of the child and of the
# reference on the arm), so a walker change cannot silently turn it into another case.
# The verdicts and causes are the library's to derive (test_traverse_stitch).
STITCH_DIR = 'stitch'


def _stitch_request(seed_id, sequence, strategy, labels=None):
    seed = {'seed_id': seed_id, 'sequence': sequence}
    if labels is not None:
        seed['labels'] = labels
    return {'seeds': [seed], 'strategy': strategy}


_KEEP_FORBID = {'direction': 'right', 'branching': {'on_reconverge': 'keep'},
                'bounds': {'max_extension_bp': 20}}
_LOOP = {'direction': 'right', 'branching': {'max_label_branches': 1, 'on_reconverge': 'keep'},
         'bounds': {'max_extension_bp': 115}, 'output': {'continuation_bp': 40}}
_WHITELIST = {'direction': 'right', 'labels': {'seed_label_kind': 'column'},
              'bounds': {'max_extension_bp': 30}}
_CHAIN = {'direction': 'right',
          'labels': {'extra': ['B', 'C', 'D'], 'change_cost': {'model': 'constant', 'value': 1},
                     'loss_budget': 3},
          'branching': {'max_label_branches': 1}, 'bounds': {'max_extension_bp': 50}}
_BRANCHY = {'direction': 'right', 'labels': {'seed_label_kind': 'column'},
            'branching': {'on_reconverge': 'keep'}, 'bounds': {'max_extension_bp': 10}}


def _with(base, **blocks):
    """|base| with each block's fields updated (a block given as a non-dict replaces it)."""
    s = copy.deepcopy(base)
    for block, values in blocks.items():
        if isinstance(values, dict):
            s.setdefault(block, {}).update(values)
        else:
            s[block] = values
    return s


STITCH_FIXTURES = {
    # ---- (a) the whitelist cases (every rule holds: qualified with exact_gated in 1a-S)
    'keep_forbid': {
        'index': 'fork',
        'doc': 'keep and forbid: the right arm diverges at 0 (R1: acc1, acc3; R2: acc2), both '
               'walks end at the radius 20; the child continues walk 0 to the end of R1. A '
               'sibling continues walk 1',
        'parent': _stitch_request('keep_forbid', ['E'], _KEEP_FORBID),
        'arm': 'right', 'walk': 0, 'bp': 40, 'sibling': 1,
        'expect': {'L': 20, 'n': 20, 'child': [20], 'sibling': [20], 'reference': [40, 40]},
    },
    'merge_only': {
        'index': 'whitelist',
        'doc': 'merge with w the arm\'s only open walk: a.fa splits ambiguously at 0 '
               '(unlimited), its X2 walk dead-ends at 10, the X1 walk (a.fa, b.fa) is open '
               'at the radius',
        'parent': _stitch_request('merge_only', ['S'], _with(
            _WHITELIST, branching={'max_label_branches': 'unlimited'}), ['a.fa', 'b.fa']),
        'arm': 'right', 'walk': None, 'bp': 70,
        'expect': {'L': 30, 'n': 30, 'child': [30], 'reference': [10, 60]},
    },
    'unlimited': {
        'index': 'whitelist',
        'doc': 'an unlimited branch allowance under keep: a.fa used one branch, b.fa none',
        'parent': _stitch_request('unlimited', ['S'], _with(
            _WHITELIST, branching={'max_label_branches': 'unlimited', 'on_reconverge': 'keep'}),
            ['a.fa', 'b.fa']),
        'arm': 'right', 'walk': None, 'bp': 70,
        'expect': {'L': 30, 'n': 30, 'child': [30], 'reference': [10, 60]},
    },
    'equal_finite': {
        'index': 'whitelist',
        'doc': 'an equal finite allowance: a.fa alone, max_label_branches 1, its one branch '
               'spent at 0; the child\'s allowance is 1 - 1 = 0',
        'parent': _stitch_request('equal_finite', ['S'], _with(
            _WHITELIST, branching={'max_label_branches': 1}), ['a.fa']),
        'arm': 'right', 'walk': None, 'bp': 70,
        'expect': {'L': 30, 'n': 30, 'child': [30], 'reference': [10, 60]},
    },
    'chain3': {
        'index': 'whitelist',
        'doc': 'a chain of three parts on one walk: b.fa alone (no split), the parent to 15, '
               'the child to 30, the grandchild to the end of b1 at 45',
        'parent': _stitch_request('chain3', ['S'], _with(
            _WHITELIST, bounds={'max_extension_bp': 15}), ['b.fa']),
        'arm': 'right', 'walk': 0, 'bp': 15, 'grandchild': {'walk': 0, 'bp': 30},
        'expect': {'L': 15, 'n': 15, 'child': [15], 'grandchild': [15], 'reference': [45]},
    },
    'fork_right': {
        'index': 'fork',
        'doc': 'the fork fixture\'s request (both arms, radius 10 < k): the continuation of '
               'right walk 0 reaches 5 bases into the seed (L 10, n 15)',
        'parent': _stitch_request('fork', ['E'], {'direction': 'both',
                                                  'bounds': {'max_extension_bp': 10}}),
        'arm': 'right', 'walk': 0, 'bp': 40,
        'expect': {'L': 10, 'n': 15, 'child': [30], 'reference': [40, 40]},
    },
    'fork_left': {
        'index': 'fork',
        'doc': 'as fork_right, on the left arm',
        'parent': _stitch_request('fork', ['E'], {'direction': 'both',
                                                  'bounds': {'max_extension_bp': 10}}),
        'arm': 'left', 'walk': 0, 'bp': 40,
        'expect': {'L': 10, 'n': 15, 'child': [30], 'reference': [40, 40]},
    },
    # ---- (b) one fixture per cause
    'sibling_merge': {
        'index': 'bubbles',
        'doc': 'merge: both alleles of the first SNP bubble end at the radius 30, each with a '
               'continuation, and reconverge at 36, beyond it (sibling_merge_unknown)',
        'parent': _stitch_request('sibling_merge', ['S'], {
            'direction': 'right', 'labels': {'seed_label_kind': 'column'},
            'branching': {'max_label_branches': 'unlimited'},
            'bounds': {'max_extension_bp': 30}}),
        'arm': 'right', 'walk': 0, 'bp': 70,
        'expect': {'L': 30, 'n': 30, 'child': [70], 'reference': [100]},
    },
    'trace': {
        'index': 'fork',
        'doc': 'keep_forbid under support trace (trace_restart)',
        'parent': _stitch_request('trace', ['E'], _with(_KEEP_FORBID, support='trace')),
        'arm': 'right', 'walk': 0, 'bp': 40,
        'expect': {'L': 20, 'n': 20, 'child': [20], 'reference': [40, 40]},
    },
    'quorum_changed': {
        'index': 'fork',
        'doc': 'keep_forbid with the child\'s quorum raised (min_successor_labels 2): '
               'strategy_differs and label_count_knob',
        'parent': _stitch_request('quorum_changed', ['E'], _KEEP_FORBID),
        'arm': 'right', 'walk': 0, 'bp': 40,
        'overrides': {'branching': {'min_successor_labels': 2}},
        'expect': {'L': 20, 'n': 20, 'child': [0], 'reference': [40, 40]},
    },
    'reset_branches': {
        'index': 'whitelist',
        'doc': 'equal_finite\'s parent with the child\'s allowance restarted '
               '(reset_branches=True): 1, not 1 - 1 (branches_unequal)',
        'parent': _stitch_request('reset_branches', ['S'], _with(
            _WHITELIST, branching={'max_label_branches': 1}), ['a.fa']),
        'arm': 'right', 'walk': None, 'bp': 70, 'reset_branches': True,
        'expect': {'L': 30, 'n': 30, 'child': [30], 'reference': [10, 60]},
    },
    'switching': {
        'index': 'chain',
        'doc': 'A -> B by a switch at 30 (constant cost 1, loss_budget 3); the walk ends at '
               'the radius 50 under B (switching)',
        'parent': _stitch_request('switching', ['S'], _CHAIN, ['A']),
        'arm': 'right', 'walk': 0, 'bp': 50,
        'expect': {'L': 50, 'n': 20, 'child': [40, 40], 'reference': [90, 90]},
    },
    'extra_foreign': {
        'index': 'chain',
        'doc': 'switching\'s parent with labels.extra B, C; the child names D, outside the '
               'parent\'s pool (extra_foreign)',
        'parent': _stitch_request('extra_foreign', ['S'], _with(
            _CHAIN, labels={'extra': ['B', 'C']}), ['A']),
        'arm': 'right', 'walk': 0, 'bp': 50,
        'overrides': {'labels': {'extra': ['A', 'C', 'D']}},
        'expect': {'L': 50, 'n': 20, 'child': [40, 40], 'reference': [90]},
    },
    'not_continued': {
        'index': 'switch_tail',
        'doc': 'z.fa is lost at 20 and switches to y.fa; at the radius 25 y.fa is alive but '
               'does not cover the continuation\'s 15 bp tail, so the child does not seed it '
               '(not_continued)',
        'parent': _stitch_request('not_continued', ['S'], {
            'direction': 'right',
            'labels': {'seed_label_kind': 'column', 'extra': ['y.fa'],
                       'change_cost': {'model': 'constant', 'value': 1}, 'loss_budget': 2},
            'bounds': {'max_extension_bp': 25}}, ['x.fa', 'z.fa']),
        'arm': 'right', 'walk': 0, 'bp': 35,
        'expect': {'L': 25, 'n': 15, 'child': [35], 'reference': [60]},
    },
    'loop_branch': {
        'index': 'loop',
        'doc': 'the external review\'s counterexample (G.4): forbid, max_label_branches 1, '
               'keep, radius 115, continuation_bp 40. The walk P U C U ends at the radius at '
               'u, its one branch spent at the first u; the child (allowance 0) sees two exits '
               'and stops at 0 bp with branch, the reference goes on to 145 (candidate_changed)',
        'parent': _stitch_request('loop_branch', ['S'], _LOOP),
        'arm': 'right', 'walk': None, 'bp': 885,
        'expect': {'L': 115, 'n': 40, 'child': [0], 'reference': [85, 145]},
    },
    'loop_branch_cut0': {
        'index': 'loop',
        'doc': 'loop_branch\'s child run with output.max_branch_events 0: its one branch '
               'event is not kept, its boundary is 0 and its walk 0 bp (decisions_unrecorded)',
        'parent': _stitch_request('loop_branch_cut0', ['S'], _LOOP),
        'arm': 'right', 'walk': None, 'bp': 885,
        'overrides': {'output': {'max_branch_events': 0}},
        'expect': {'L': 115, 'n': 40, 'child': [0], 'reference': [85, 145]},
    },
    'loop_same_stop': {
        'index': 'loop',
        'doc': 'loop_branch with the default continuation_bp: the child is seeded with the '
               'whole walk, so it blocks the loop\'s first edge as rejoined_seed, an edge the '
               'parent used (seed_reentry, same_stop), and goes on along Y as the reference does',
        'parent': _stitch_request('loop_same_stop', ['S'], {
            'direction': 'right',
            'branching': {'max_label_branches': 1, 'on_reconverge': 'keep'},
            'bounds': {'max_extension_bp': 115}}),
        'arm': 'right', 'walk': None, 'bp': 885,
        'expect': {'L': 115, 'n': 115, 'child': [30], 'reference': [85, 145]},
    },
    'loop_interior': {
        'index': 'loop',
        'doc': 'unlimited branches, radius 90, continuation_bp 15: the walk into C ends inside '
               'it; 25 bp into the child it reaches u, where the loop\'s first edge, used by '
               'the parent, is admitted again: the walk through it takes a step in B '
               '(rule_violation), the walk through Y has a changed candidate '
               '(candidate_changed), both past the child\'s first base',
        'parent': _stitch_request('loop_interior', ['S'], _with(
            _LOOP, branching={'max_label_branches': 'unlimited', 'on_reconverge': 'keep'},
            bounds={'max_extension_bp': 90}, output={'continuation_bp': 15})),
        'arm': 'right', 'walk': None, 'bp': 910,
        'expect': {'L': 90, 'n': 15, 'child': [55, 59], 'reference': [85, 145]},
    },
    'earlier_stop': {
        'index': 'repeat',
        'doc': 'unlimited branches, keep, radius 70: the walk along Q is the child\'s seed; the '
               'child blocks the new edge into R\'s first k-mer as rejoined_seed, where the '
               'reference takes it and stops one base later with edge_reuse (seed_reentry, '
               'earlier_stop)',
        'parent': _stitch_request('earlier_stop', ['S'], {
            'direction': 'right',
            'branching': {'max_label_branches': 'unlimited', 'on_reconverge': 'keep'},
            'bounds': {'max_extension_bp': 70}}),
        'arm': 'right', 'walk': None, 'bp': 930,
        'expect': {'L': 70, 'n': 70, 'child': [30], 'reference': [66, 101]},
    },
    'splits45': {
        'index': 'splits',
        'doc': 'the full Fable review\'s split-count case (G.5): k 11, forbid, merge, '
               'max_label_branches 0, max_splits_per_path 2, radius 45; the walk\'s leaf holds '
               'two splits through a merge parent, one on its displayed chain; the child '
               'splits again, which the reference refuses (split_limit_restarted)',
        'parent': _stitch_request('splits45', ['S'], {
            'direction': 'right', 'labels': {'seed_label_kind': 'column'},
            'branching': {'max_label_branches': 0, 'max_splits_per_path': 2},
            'bounds': {'max_extension_bp': 45}}),
        'arm': 'right', 'walk': None, 'bp': 40,
        'expect': {'L': 45, 'n': 14, 'child': [40, 40], 'reference': [18, 50]},
    },
    'events_cut': {
        'index': 'branchy',
        'doc': 'max_label_branches 2: the child splits ambiguously at 30 and again at 50, and '
               'keeps one branch event (output.max_branch_events 1): its boundary is 50 and '
               'its walks end beyond it (decisions_unrecorded)',
        'parent': _stitch_request('events_cut', ['S'], _with(
            _BRANCHY, branching={'max_label_branches': 2})),
        'arm': 'right', 'walk': 0, 'bp': 190,
        'overrides': {'output': {'max_branch_events': 1}},
        'expect': {'L': 10, 'n': 15, 'child': [40, 80, 90], 'reference': [50, 90, 100]},
    },
    'events_at_end': {
        'index': 'branchy',
        'doc': 'max_label_branches 1: the child splits at 30 and its label ends with branch at '
               '50, where the one branch event kept is followed by one not kept: the boundary '
               '50 equals the walk\'s end (decisions_unrecorded)',
        'parent': _stitch_request('events_at_end', ['S'], _with(
            _BRANCHY, branching={'max_label_branches': 1})),
        'arm': 'right', 'walk': 0, 'bp': 190,
        'overrides': {'output': {'max_branch_events': 1}},
        'expect': {'L': 10, 'n': 15, 'child': [50, 40], 'reference': [60, 50]},
    },
    'circular': {
        'index': 'cycle',
        'doc': 'the GPT review\'s circular record (G.6): k 15, radius 25, continuation_bp 15; '
               'the parent\'s continuation is the child\'s seed and the child\'s continuation '
               'the parent\'s, and the child\'s walk enters the root seed\'s node '
               '(rule_violation)',
        'parent': _stitch_request('circular', [('=', 'GGCCCCGCCGCAACT')], {
            'direction': 'right', 'branching': {'on_reconverge': 'keep'},
            'bounds': {'max_extension_bp': 25}, 'output': {'continuation_bp': 15}}),
        'arm': 'right', 'walk': 0, 'bp': 25,
        'expect': {'L': 25, 'n': 15, 'child': [25], 'reference': [49]},
    },
    'no_bases': {
        'index': 'fork',
        'doc': 'keep_forbid with output.sequences false: the rule check cannot run (no_bases)',
        'parent': _stitch_request('no_bases', ['E'], _with(_KEEP_FORBID,
                                                          output={'sequences': False})),
        'arm': 'right', 'walk': 0, 'bp': 40,
        'expect': {'L': 20, 'n': 20, 'child': [20], 'reference': [40, 40]},
    },
}

# CLI fixtures of the request budgets (DESIGN-traverse-graphlet.md §14), written to
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
    them, plus output.detail/timing and clamped."""
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


def _merge_overrides(strategy, overrides):
    """next_request()'s deep merge of keyword overrides (a change_cost naming another model
    replaces the whole; the fixtures override nothing else that needs care)."""
    for k, v in overrides.items():
        cur = strategy.get(k)
        if isinstance(v, dict) and isinstance(cur, dict) and not (
                k == 'change_cost' and v.get('model', cur.get('model')) != cur.get('model')):
            _merge_overrides(cur, v)
        else:
            strategy[k] = copy.deepcopy(v)


def stitch_child_request(parent, full, arm, walk, bp, overrides=None, reset_branches=False):
    """The request parent.next_request(arm, [walk], bp, reset_branches=..., **overrides)
    builds (ops.next_request), for the cases the fixtures use: one walk, the parent's
    normalized strategy with direction = |arm| and max_extension_bp = |bp|, the walk's
    continuation as the seed with its labels named (constrain), labels.loss_budget reduced
    by the largest terminal loss of those labels, branching.max_label_branches by their
    largest terminal branch count (unless "unlimited", 0, or reset_branches), labels.extra
    rebuilt around the new seed labels (forbid: none; constant: every other label of the
    retrieval when one switch fits the budget), then the overrides merged. |parent| is the
    detail graphlet response, |full| the detail full one of the same request. A cost table
    is not supported here (the library's switch search would be needed)."""
    overrides = copy.deepcopy(overrides or {})
    strategy = copy.deepcopy(parent['strategy'])
    strategy.pop('clamped', None)
    strategy['direction'] = arm
    strategy.setdefault('bounds', {})['max_extension_bp'] = bp
    r = full['results'][0]
    names = [l['name'] for l in r['label_dict']]
    path = r['arms'][arm]['paths'][walk]
    cont = path['continuation']
    seed = {'sequence': cont['sequence']}
    labels = strategy['labels']
    if labels['mode'] == 'constrain':
        if cont['labels']:
            seed['labels'] = [names[l] for l in cont['labels']]
        ends = {e['label']: e for e in path.get('end_labels', [])}
        budget = float(labels.get('loss_budget') or 0.0)
        top = max([ends[l]['loss'] for l in cont['labels'] if l in ends] or [0.0])
        effective = budget
        if budget > 0 and top > 0:
            effective = max(0.0, budget - top)
            labels['loss_budget'] = effective
        cost = labels.get('change_cost') or {'model': 'forbid'}
        given = overrides.get('labels') or {}
        if 'change_cost' in given:
            cost = dict(cost, **given['change_cost']) \
                if given['change_cost'].get('model', cost['model']) == cost['model'] \
                else given['change_cost']
        model = cost.get('model', 'forbid')
        if model == 'table':
            raise NotImplementedError('a cost table: the fixtures use forbid or constant')
        seeded = set(cont['labels'])
        if model == 'constant' and cost.get('value', 0) <= given.get('loss_budget', effective):
            labels['extra'] = [n for i, n in enumerate(names) if i not in seeded]
        else:
            labels['extra'] = []
        branches = max([cont.get('branches_used', 0)]
                       + [ends[l]['branches'] for l in cont['labels'] if l in ends])
        original = strategy['branching'].get('max_label_branches')
        if original == 'unlimited' or (isinstance(original, int) and original > 0):
            strategy['branching']['max_label_branches'] = \
                original if original == 'unlimited' or reset_branches \
                else max(0, original - branches)
    _merge_overrides(strategy, overrides)
    out = strategy.get('output') or {}
    if out.get('coordinates') is not True:
        out.pop('max_coordinate_occurrences', None)
    request = {'seeds': [seed]}
    if parent.get('release'):
        request['release'] = parent['release']
    request['strategy'] = strategy
    return request


def _arm_lengths(resp, arm):
    a = resp['results'][0]['arms'].get(arm)
    return None if a is None else [p['length_bp'] for p in a['paths']]


def _pick_walk(full, arm, walk):
    """|walk| as given, or (None) the first walk of the arm that ended at the radius."""
    if walk is not None:
        return walk
    for p in full['results'][0]['arms'][arm]['paths']:
        if p.get('path_reason') == 'max_extension_bp':
            return p['id']
    raise RuntimeError('no walk of the %s arm ended at the radius' % arm)


def build_stitch(metagraph, root, name):
    """-> ({filename: text} for stitching fixture |name| (in STITCH_DIR), [what differs from
    its expect])."""
    fx = STITCH_FIXTURES[name]
    spec = INDEXES[fx['index']]
    graph = os.path.join(root, fx['index'], 'graph.dbg')
    if not os.path.exists(graph):
        build_index(fx['index'], metagraph, root)
    annotation = annotation_of(fx['index'], root)
    arm = fx['arm']

    def run(request):
        """-> (detail graphlet response, detail full response) of |request|."""
        outs = {}
        for detail in ('graphlet', 'full'):
            r = copy.deepcopy(request)
            r['strategy'].setdefault('output', {})['detail'] = detail
            outs[detail] = run_traverse(metagraph, graph, annotation, r)
        return outs['graphlet'], outs['full']

    parent_req = materialize(fx['parent'], blocks_of(spec))
    parent, parent_full = run(parent_req)
    walk = _pick_walk(parent_full, arm, fx['walk'])
    child_req = stitch_child_request(parent, parent_full, arm, walk, fx['bp'],
                                     fx.get('overrides'), fx.get('reset_branches', False))
    child, child_full = run(child_req)
    reference_req = copy.deepcopy(parent_req)
    radius = parent['strategy']['bounds']['max_extension_bp']
    reference_req['strategy'].setdefault('bounds', {})['max_extension_bp'] = radius + fx['bp']
    reference, reference_full = run(reference_req)
    requests = {'fixture_index': {'cli': fx['index'], 'k': spec['k'], 'doc': fx['doc'],
                                  'arm': arm, 'walk': walk},
                'parent': parent_req, 'child': child_req, 'reference': reference_req}
    responses = {'parent': parent, 'child': child, 'reference': reference}
    fulls = {'parent': parent_full, 'child': child_full, 'reference': reference_full}
    if 'sibling' in fx:
        sib_req = stitch_child_request(parent, parent_full, arm, fx['sibling'], fx['bp'],
                                       fx.get('overrides'), fx.get('reset_branches', False))
        requests['fixture_index']['sibling'] = fx['sibling']
        requests['sibling'] = sib_req
        responses['sibling'], fulls['sibling'] = run(sib_req)
    if 'grandchild' in fx:
        gc = fx['grandchild']
        gc_req = stitch_child_request(child, child_full, arm, gc['walk'], gc['bp'])
        requests['fixture_index']['grandchild'] = gc['walk']
        requests['grandchild'] = gc_req
        responses['grandchild'], fulls['grandchild'] = run(gc_req)
        reference_req['strategy']['bounds']['max_extension_bp'] += gc['bp']
        responses['reference'], fulls['reference'] = run(reference_req)
    have = {'L': parent_full['results'][0]['arms'][arm]['paths'][walk]['length_bp'],
            'n': len(child_req['seeds'][0]['sequence'])}
    for role in ('child', 'sibling', 'grandchild', 'reference'):
        if role in fulls:
            have[role] = _arm_lengths(fulls[role], arm)
    problems = ['%s: %s is %r, expected %r' % (name, k, have.get(k), v)
                for k, v in fx['expect'].items() if have.get(k) != v]
    files = {name + '.request.json': dumps(requests)}
    for role, resp in responses.items():
        files['%s.%s.json' % (name, role)] = dumps(resp)
    return {os.path.join(STITCH_DIR, f): text for f, text in files.items()}, problems, have


# The CLI fixtures are the reference platform's answers: Linux, where the servers run (the
# CI's Linux build writes and checks them byte for byte). The memory account counts the
# platform's object sizes (libstdc++ and libc++ lay out std::string and the structs holding
# one differently; SPEC-labeled-traversal-core.md section 6.8), so another platform states other
# byte counts and stops a memory-bounded walk elsewhere. There --check compares every byte but
# the ACCOUNTING_FIELDS, and checks a memory-stopped answer by what it states
# (memory_stop_problems) instead of its bytes.
REFERENCE_PLATFORM = 'Linux'
# the modelled byte counts of a document: (path from the document, key), '[]' every element
ACCOUNTING_FIELDS = (
    (('work',), 'memory_bytes'),
    (('results', '[]', 'work'), 'memory_bytes'),
)
CHECK_MODE = None   # 'reference' (byte for byte) or 'portable', set by main()


def without_accounting(doc):
    """|doc| without its ACCOUNTING_FIELDS"""
    doc = copy.deepcopy(doc)

    def drop(node, path, key):
        if not path:
            if isinstance(node, dict):
                node.pop(key, None)
            return
        head, rest = path[0], path[1:]
        if head == '[]':
            for x in node if isinstance(node, list) else []:
                drop(x, rest, key)
        elif isinstance(node, dict) and head in node:
            drop(node[head], rest, key)
    for path, key in ACCOUNTING_FIELDS:
        drop(doc, path, key)
    return doc


def memory_stopped(doc):
    """whether a document's walk was stopped by its memory budget"""
    return any((r.get('resource_stop') or {}).get('resource') == 'memory'
               for r in doc.get('results', []))


def memory_stop_problems(have, new):
    """A memory-stopped answer of a platform other than the reference, checked by what it
    states: the same stop (resource, phase, scope, budget, actions) and outcome as the
    reference's; each arm's complete_to_bp where the stop's message says the walk stopped and
    where its cap_trigger stands, no longer than its walks; its output_bp at least that and
    equal to its bases; its peak within the budget plus the soft excess it states. -> the
    statements that do not hold (empty: the answer is true to itself)"""
    out = []
    if len(have.get('results', [])) != len(new.get('results', [])):
        return ['a different number of results']
    for i, (a, b) in enumerate(zip(have['results'], new['results'])):
        sa, sb = a.get('resource_stop') or {}, b.get('resource_stop') or {}
        for k in ('resource', 'phase', 'scope', 'requested', 'actions'):
            if sa.get(k) != sb.get(k):
                out.append(f'results[{i}].resource_stop.{k}: {sb.get(k)!r}, '
                           f'the reference states {sa.get(k)!r}')
        if a.get('outcome') != b.get('outcome'):
            out.append(f'results[{i}].outcome: {b.get("outcome")!r}')
        m = re.search(r'stopped the walk at (\d+) bp on the (\w+) arm', sb.get('message', ''))
        for name, arm in (b.get('arms') or {}).items():
            done = arm.get('complete_to_bp')
            trigger = arm.get('cap_trigger') or {}
            counts, counters = arm.get('counts') or {}, arm.get('counters') or {}
            if trigger and trigger.get('at_bp') != done:
                out.append(f'{name} arm: cap_trigger.at_bp {trigger.get("at_bp")} '
                           f'!= complete_to_bp {done}')
            if m and m.group(2) == name and int(m.group(1)) != done:
                out.append(f'{name} arm: the stop says {m.group(1)} bp, complete_to_bp {done}')
            if counts.get('max_bp') is not None and done is not None and done > counts['max_bp']:
                out.append(f'{name} arm: complete_to_bp {done} > max_bp {counts["max_bp"]}')
            if 'output_bp' in counters and done is not None and counters['output_bp'] < done:
                out.append(f'{name} arm: output_bp {counters["output_bp"]} < complete_to_bp')
            if 'bases' in counts and 'output_bp' in counters \
                    and counts['bases'] != counters['output_bp']:
                out.append(f'{name} arm: bases {counts["bases"]} != output_bp')
        budget = sb.get('requested')
        soft = max([l.get('observed') or 0 for l in b.get('limitations', [])
                    if l.get('kind') == 'memory_bound_soft'] or [0])
        peak = (b.get('work') or {}).get('memory_bytes')
        if budget is not None and peak is not None and peak > (budget + soft) * 2**20:
            out.append(f'results[{i}]: peak {peak} bytes over the budget {budget} MiB '
                       f'and the soft excess {soft} MiB')
    return out


def mgt_statements(text):
    """a graphlet's records that state the walk's frame and stop: its header, seed, labels
    and options lines whole, the stop (Q) by its scope, resource and phase, the limitations
    (K) by their kind"""
    keep = []
    for line in text.split('\n'):
        fields = line.split(' ')
        if fields[0] in ('H', 'S', 'L', 'O'):
            keep.append(line)
        elif fields[0] == 'Q':
            keep.append(' '.join(fields[:4]))
        elif fields[0] == 'K':
            keep.append(' '.join(fields[:3]))
    return sorted(keep)


def differs(path, text):
    if not os.path.exists(path):
        return True
    have = read_text(path)
    if have == text:
        return False
    if path.endswith('.json'):
        old, new = json.loads(have), json.loads(text)
        if CHECK_MODE == 'portable':
            if memory_stopped(old):
                problems = memory_stop_problems(old, new)
                for p in problems:
                    print(f'{os.path.basename(path)}: {p}', file=sys.stderr)
                return bool(problems)
            return without_accounting(old) != without_accounting(new)
        return old != new
    if CHECK_MODE == 'portable' and path.endswith('.mgt'):
        companion = path[:-len('.mgt')] + '.graphlet.json'
        if os.path.exists(companion) and memory_stopped(json.loads(read_text(companion))):
            return mgt_statements(have) != mgt_statements(text)
    return True


def show_diff(path, text, limit=40):
    """--check: the first |limit| lines of a unified diff of the committed fixture |path|
    against the regenerated |text| (a JSON document both re-indented, its keys in their
    order), on stderr, so that a failed check shows what differs"""
    have = read_text(path) if os.path.exists(path) else ''
    a, b = have, text
    if path.endswith('.json') and have:
        try:
            a = json.dumps(json.loads(have), indent=1, ensure_ascii=False) + '\n'
            b = json.dumps(json.loads(text), indent=1, ensure_ascii=False) + '\n'
        except ValueError:
            pass
    lines = list(difflib.unified_diff(a.splitlines(keepends=True), b.splitlines(keepends=True),
                                      'committed/' + os.path.basename(path),
                                      'regenerated/' + os.path.basename(path), n=2))
    sys.stderr.writelines(lines[:limit])
    if len(lines) > limit:
        print('... (%d more diff lines)' % (len(lines) - limit), file=sys.stderr)
    elif lines and not lines[-1].endswith('\n'):
        print(file=sys.stderr)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--dir', default=DEFAULT_DIR)
    ap.add_argument('--check', action='store_true', help='verify, write nothing')
    ap.add_argument('--from-cli', metavar='METAGRAPH',
                    help='regenerate (or with --check verify) the CLI fixtures')
    ap.add_argument('--only', nargs='*')
    ap.add_argument('--keep', metavar='DIR', help='build the indexes here and keep them')
    ap.add_argument('--platform-mode', choices=('auto', 'reference', 'portable'),
                    default='auto',
                    help='--check: byte for byte (reference: the Linux answers the fixtures are) '
                         'or portable (every byte but the accounting fields, a memory-stopped '
                         'answer by its statements); auto: reference on Linux')
    args = ap.parse_args(argv)
    global CHECK_MODE
    CHECK_MODE = args.platform_mode
    if CHECK_MODE == 'auto':
        CHECK_MODE = 'reference' if platform.system() == REFERENCE_PLATFORM else 'portable'
    if not args.from_cli:
        # the hand-made documents are written by this script, the same on every platform
        CHECK_MODE = 'reference'
    if args.check and args.from_cli:
        if CHECK_MODE == 'reference':
            print('check mode: reference (byte for byte, the %s answers)' % REFERENCE_PLATFORM,
                  flush=True)
        else:
            print('check mode: portable on %s (byte for byte but %s; a memory-stopped answer '
                  'by its statements)' % (platform.system(), ', '.join(
                      '.'.join(path + (key,)) for path, key in ACCOUNTING_FIELDS)), flush=True)
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
                        show_diff(path, text)
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
                        show_diff(path, text)
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
                    show_diff(path, text)
                else:
                    with open(path, 'w', encoding='utf-8', newline='') as f:
                        f.write(text)
                    print('wrote', fname)
            for name in STITCH_FIXTURES:
                if args.only and name not in args.only:
                    continue
                files, problems, _ = build_stitch(args.from_cli, root, name)
                if problems:
                    stale.extend(problems)
                    continue
                os.makedirs(os.path.join(args.dir, STITCH_DIR), exist_ok=True)
                for fname, text in files.items():
                    path = os.path.join(args.dir, fname)
                    if not differs(path, text):
                        continue
                    if args.check:
                        stale.append(fname)
                        show_diff(path, text)
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
                    show_diff(path, text)
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
