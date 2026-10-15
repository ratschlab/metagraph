"""Build the fixture indexes with the engine's own CLI (ORACLE-DESIGN.md section 2).

The pipeline is the integration tests' and build_mini_refseq.sh's: `build` (with the dummy-edge
mask), `annotate --anno-filename` (with `--coordinates` on BASIC), on BASIC the row_diff stages
0/1/2, `row_diff_brwt_coord --greedy` and `relax_brwt`, the columns renamed to the fixture's
column names, and on BASIC the CoordToHeader sidecar (`.seqs`) in the final column order. Every
command and its output go to <out>/build/logs. Nothing is ever deleted.
"""

import os
import subprocess

from . import fixtures

# config -> (graph mode, k, masked, annotation kind, fixture set)
CONFIGS = {
    'basic_k11': ('basic', 11, True, 'row_diff_brwt_coord', 'records'),
    'basic_k11_unmasked': ('basic', 11, False, 'row_diff_brwt_coord', 'records'),
    'basic_k11_unchecked': ('basic', 11, False, 'row_diff_brwt_coord', 'records'),
    'canonical_k11': ('canonical', 11, True, 'column', 'records'),
    'primary_k11': ('primary', 11, True, 'column', 'records'),
    'basic_k15': ('basic', 15, True, 'row_diff_brwt_coord', 'records'),
    # even k: k-mers equal to their reverse complement (one row on CANONICAL and PRIMARY; the
    # palindrome scans of an even-k wrapped PRIMARY graph, SPEC 7.4-7.6)
    'basic_k10': ('basic', 10, True, 'row_diff_brwt_coord', 'records'),
    'canonical_k10': ('canonical', 10, True, 'column', 'records'),
    'primary_k10': ('primary', 10, True, 'column', 'records'),
    # phase 2: the traversal fixture (traverse_fixtures.py), for /traverse only
    'trav_basic_k11': ('basic', 11, True, 'row_diff_brwt_coord', 'trav'),
    'trav_canonical_k11': ('canonical', 11, True, 'column', 'trav'),
    'trav_primary_k11': ('primary', 11, True, 'column', 'trav'),
}
# the server flags every config gets: the floor lowered (ORACLE-DESIGN.md section 2) and ample
# caps, so that no generated request is cut by a server default unless it asks to be
SERVER_FLAGS = ['--pattern-min-information-bits', '0', '--pattern-max-contexts', '100000',
                '--pattern-max-occurrences', '100000']
EXTRA_FLAGS = {'basic_k11_unchecked': ['--pattern-max-checked-entries', '0']}


class BuildError(RuntimeError):
    pass


class Builder:
    def __init__(self, metagraph, out_dir):
        self.mg = metagraph
        self.root = os.path.join(out_dir, 'build')
        self.logs = os.path.join(self.root, 'logs')
        os.makedirs(self.logs, exist_ok=True)
        self.n = 0

    def run(self, args, cwd=None, name='step'):
        self.n += 1
        log = os.path.join(self.logs, f'{self.n:02d}_{name}.log')
        with open(log, 'w') as f:
            f.write(' '.join([self.mg] + args) + '\n\n')
            f.flush()
            res = subprocess.run([self.mg] + args, cwd=cwd, stdout=f, stderr=subprocess.STDOUT,
                                 stdin=subprocess.DEVNULL)
        if res.returncode != 0:
            raise BuildError(f'{name} failed (exit {res.returncode}); see {log}')
        return log

    def output(self, args, cwd=None):
        res = subprocess.run([self.mg] + args, cwd=cwd, stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, text=True)
        if res.returncode != 0:
            raise BuildError(f'{" ".join(args)} failed: {res.stdout[-2000:]}')
        return res.stdout


def column_names(builder, annotation, graph):
    """The column names of |annotation| in column order (`stats --print-col-names`)."""
    text = builder.output(['stats', '-a', annotation, '--print-col-names', graph])
    lines = text.splitlines()
    for i, line in enumerate(lines):
        if line.startswith('Number of columns:'):
            n = int(line.split(':')[1])
            return [l.strip() for l in lines[i + 1:i + 1 + n]]
    raise BuildError(f'no column list in stats of {annotation}')


def _rename(builder, d, annotation, graph, out_base, name):
    rules = os.path.join(d, 'rename_cols.txt')
    with open(rules, 'w') as f:
        for col in column_names(builder, annotation, graph):
            f.write(f'{col} {os.path.splitext(os.path.basename(col))[0]}\n')
    builder.run(['transform_anno', '--rename-cols', rules, '-o', out_base, annotation], cwd=d,
                name=name)


def build_fasta(out_dir, columns, name='fasta'):
    d = os.path.join(out_dir, 'build', name)
    os.makedirs(d, exist_ok=True)
    return fixtures.write_fasta(columns, d)


def build_basic(builder, fasta, k, mode='basic', prefix='', order=None):
    """A BASIC (or CANONICAL with |mode|) masked graph with its annotation; -> paths. |order|:
    the columns in annotation order (their column ids)."""
    d = os.path.join(builder.root, f'{prefix}{mode}_k{k}')
    os.makedirs(d, exist_ok=True)
    files = [fasta[c] for c in (order or fixtures.COLUMN_ORDER)]
    graph_base = os.path.join(d, 'graph')
    graph = graph_base + '.dbg'
    builder.run(['build', '-p', '1', '--mode', mode, '--graph', 'succinct', '-k', str(k),
                 '--mask-dummy', '--in-ram', '-o', graph_base] + files, name=f'{mode}{k}_build')
    anno_base = os.path.join(d, 'annotation')
    if mode != 'basic':
        builder.run(['annotate', '-p', '1', '-i', graph, '--anno-filename', '-o', anno_base]
                    + files, name=f'{mode}{k}_annotate')
        final = os.path.join(d, 'final')
        _rename(builder, d, anno_base + '.column.annodbg', graph, final, f'{mode}{k}_rename')
        return {'graph': graph, 'annotation': final + '.column.annodbg', 'dir': d}
    builder.run(['annotate', '-p', '1', '-i', graph, '--anno-filename', '--coordinates', '-o',
                 anno_base] + files, name=f'basic{k}_annotate')
    rd = os.path.join(d, 'rd')
    os.makedirs(rd, exist_ok=True)
    for stage in ('0', '1', '2'):
        builder.run(['transform_anno', '-p', '1', '--anno-type', 'row_diff', '--row-diff-stage',
                     stage, '--coordinates', '-i', graph, '-o', os.path.join(rd, 'out'),
                     anno_base + '.column.annodbg'], name=f'basic{k}_row_diff_{stage}')
    builder.run(['transform_anno', '-p', '1', '--anno-type', 'row_diff_brwt_coord', '--greedy',
                 '-i', graph, '-o', anno_base, os.path.join(rd, 'annotation.column.annodbg')],
                name=f'basic{k}_brwt_coord')
    builder.run(['relax_brwt', '-p', '1', '--relax-arity', '32', '-o', anno_base + '.relaxed',
                 anno_base + '.row_diff_brwt_coord.annodbg'], name=f'basic{k}_relax')
    final = os.path.join(d, 'final')
    _rename(builder, d, anno_base + '.relaxed.row_diff_brwt_coord.annodbg', graph, final,
            f'basic{k}_rename')
    annotation = final + '.row_diff_brwt_coord.annodbg'
    order = column_names(builder, annotation, graph)
    ordered = [fasta[c] for c in order]
    builder.run(['annotate', '-p', '1', '-i', graph, '--anno-filename', '--index-header-coords',
                 '-o', final] + ordered, name=f'basic{k}_seqs')
    if not os.path.exists(final + '.seqs'):
        raise BuildError('the .seqs sidecar was not written')
    return {'graph': graph, 'annotation': annotation, 'dir': d, 'columns': order}


def build_primary(builder, fasta, k, prefix='', order=None):
    d = os.path.join(builder.root, f'{prefix}primary_k{k}')
    os.makedirs(d, exist_ok=True)
    files = [fasta[c] for c in (order or fixtures.COLUMN_ORDER)]
    canon = os.path.join(d, 'canonical')
    builder.run(['build', '-p', '1', '--mode', 'canonical', '--graph', 'succinct', '-k', str(k),
                 '-o', canon] + files, name=f'primary{k}_canonical')
    contigs = os.path.join(d, 'primary_contigs')
    builder.run(['transform', '-p', '1', '--to-fasta', '--primary-kmers', '-o', contigs,
                 canon + '.dbg'], name=f'primary{k}_contigs')
    graph_base = os.path.join(d, 'graph')
    graph = graph_base + '.dbg'
    builder.run(['build', '-p', '1', '--mode', 'primary', '--graph', 'succinct', '-k', str(k),
                 '--mask-dummy', '--in-ram', '-o', graph_base, contigs + '.fasta.gz'],
                name=f'primary{k}_build')
    anno_base = os.path.join(d, 'annotation')
    builder.run(['annotate', '-p', '1', '-i', graph, '--anno-filename', '-o', anno_base] + files,
                name=f'primary{k}_annotate')
    final = os.path.join(d, 'final')
    _rename(builder, d, anno_base + '.column.annodbg', graph, final, f'primary{k}_rename')
    return {'graph': graph, 'annotation': final + '.column.annodbg', 'dir': d}


def unmasked_view(builder, built, name):
    """A directory of symlinks to a masked BASIC graph's files without its .edgemask: the loader
    reads a mask beside the path it is given, so the same graph is served without it."""
    d = os.path.join(builder.root, name)
    os.makedirs(d, exist_ok=True)
    src = built['dir']
    for f in os.listdir(src):
        if f.startswith('graph.dbg'):
            link = os.path.join(d, f)
            if not os.path.lexists(link):
                os.symlink(os.path.join(src, f), link)
    if os.path.exists(os.path.join(d, 'graph.edgemask')):
        raise BuildError('the unmasked view holds a mask')
    return {'graph': os.path.join(d, 'graph.dbg'), 'annotation': built['annotation'], 'dir': d}


def build_all(metagraph, out_dir, fixture_sets, configs):
    """Build what |configs| need from |fixture_sets| ({set: {column: [Record]}} in column order);
    -> {config: {'graph', 'annotation', 'server_args'}}."""
    builder = Builder(metagraph, out_dir)
    fastas = {}
    built = {}
    cache = {}
    for name in configs:
        mode, k, masked, _, fset = CONFIGS[name]
        columns = fixture_sets[fset]
        if fset not in fastas:
            fastas[fset] = build_fasta(out_dir, columns,
                                       'fasta' if fset == 'records' else f'fasta_{fset}')
        fasta = fastas[fset]
        order = list(columns)
        prefix = '' if fset == 'records' else f'{fset}_'
        if mode == 'primary':
            b = build_primary(builder, fasta, k, prefix, order)
        else:
            key = (fset, mode, k)
            if key not in cache:
                cache[key] = build_basic(builder, fasta, k, mode, prefix, order)
            b = cache[key]
            if not masked:
                b = unmasked_view(builder, b, name)
        args = ['-i', b['graph'], '-a', b['annotation']] + SERVER_FLAGS + EXTRA_FLAGS.get(name, [])
        built[name] = dict(b, server_args=args)
    return built
