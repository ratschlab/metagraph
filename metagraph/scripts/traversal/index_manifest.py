#!/usr/bin/env python3
"""Write or verify the manifest of a MetaGraph index bundle (DESIGN-traverse-graphlet.md §3.1).

The identity of an index is the content of the files a server loads -- the graph, the
annotation and their sidecars -- not its metadata: two annotations over one graph can agree on
k, columns and counts and still disagree on which column carries which sequence. The manifest
lists every file of the bundle with its size and sha256, and

    index_fp = sha256 over the lines "<path>\\t<size>\\t<sha256>\\n", sorted by path (bytes)

which is exactly what `metagraph server_query|traverse --index-manifest FILE` recomputes. The
server does not re-hash the bundle (hundreds of GB) at start-up: it checks that every file it
loads is listed with its size and that the stated index_fp is the digest of the list.

    index_manifest.py -i GRAPH -a ANNOTATION [-o MANIFEST] [--name NAME] [--extra FILE ...]
                      [--no-coord-mapping]
    index_manifest.py -i GRAPH -a ANNOTATION --inventory   # the files the server loads (JSON:
                                                 # files = the identity, derived = mask, Bloom)
    index_manifest.py --verify MANIFEST          # re-hash every listed file and compare
    index_manifest.py --server-csv GRAPHS.csv [--jobs N] [--digests FILE ...] [--digests-only]
                      [--match-base-names] [--out-dir DIR] [--write-csv OUT.csv] [--force]
                      [--no-coord-mapping]

Batch mode (--server-csv) writes one manifest per (graph, annotation) pair of a multi-graph
server's list, read as the server reads it (name,graph_path,annotation_path[,manifest_path
[,index_ns]], one line each, empty lines skipped): every distinct file is hashed once, by N
parallel streams (hashlib releases the GIL), or taken from precomputed digests (--digests: a
file of "<sha256>  <path>" or "<sha256> *<path>" lines, as sha256sum writes them, e.g. from
storage checksums; '#' comments). A digest line names a bundle file when its path, relative
paths taken from the digest file's own directory, is that file (by real path); two different
digests for one file are an error, however the paths are spelled. A digest file whose paths
are not this machine's (storage checksums) can be matched by base name with
--match-base-names: only lines whose paths name no file here, only a base name that one such
line carries (two with different digests match nothing), never one digest line for two
different files (an error: bundles often share base names such as graph.dbg, and the server,
which checks sizes only, could not notice a digest of another bundle's file), and every such
match is listed. --digests-only
refuses to hash a file without a digest. Each file entry records how its digest was obtained
("digest": "hashed" | "precomputed" | "precomputed-by-base-name", metadata the server
ignores). A pair's manifest goes to the manifest_path its lines give (two different ones are
an error), else to DIR/<name>.<annotation base>.manifest.json with --out-dir (name: the
pair's first line's), else next to the annotation as in single mode; an existing file is
kept unless --force. --write-csv writes the list back with the manifest column of every line
filled with its pair's manifest (the index_ns column kept), ready for server_query.

The bundle is the loader dependency inventory of the pair (load_inventory, which mirrors
metagraph's index_load_inventory; `metagraph traverse --index-inventory -i GRAPH -a ANNOTATION`
prints the binary's, and an integration test compares the two): the graph and the annotation,
and every identity sidecar the server loads for them, its path derived from the LISTED spelling
as the loaders derive it (next to a symlink, not next to its target) -- <graph>.anchors and
<graph>.rd_succ for a .row_diff.annodbg (required), and <annotation without
.<type>.annodbg>.seqs for a coordinate annotation unless --no-coord-mapping. The graph's DERIVED
data is not part of it (the owner's decision #17 of 2026-10-08; derived_files mirrors
metagraph's index_derived_files): <graph without .dbg>.edgemask, loaded when it opens, and
<graph without .dbg>.bloom, loaded when it exists and the mask was loaded, are computed from the
graph alone and decide which counts are exact and how fast k-mers are looked up, never what an
exact answer is; adding, removing or rebuilding them leaves index_fp unchanged. A manifest never
lists them: this tool leaves them out, refuses them with --extra, --verify flags a manifest that
lists one (any entry named *.edgemask or *.bloom), and the server refuses such a manifest.
--inventory prints both, the identity files and the derived ones (with whether the loader reads
them). Files the server does not open are not part of the bundle (a column annotation's .coords,
the graph's .weights): add them with --extra if wanted, except a further graph (*dbg) or
annotation (*.annodbg): a manifest describes one graph with one annotation, and the server
refuses one that lists another -- or an optional file of the pair's inventory that it does not
load (a .seqs missing beside the listed spelling or under --no-coord-mapping), so that index_fp
identifies the loaded files. The server refuses a manifest that does not cover a file of the
inventory. Paths in the manifest are the files' BASE NAMES, so the identity does not
depend on where the bundle or the manifest lives (the server matches loaded files by base name
too); the base names in one manifest must be distinct (the server, and --verify, refuse one
that lists a base name twice). --verify looks for the files next to the manifest, or under
--root. Stdlib only.
"""
import argparse
import hashlib
import json
import os
import re
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor



def sha256_file(path, progress=None):
    h = hashlib.sha256()
    done = 0
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 24), b''):
            h.update(block)
            done += len(block)
            if progress:
                progress(done)
    return h.hexdigest()


def annotation_base(path):
    # annotation.relaxed.row_diff_brwt.annodbg -> annotation.relaxed (the CoordToHeader
    # sidecar is <base>.seqs, written next to the annotation by `annotate --index-header-coords`)
    name = path[:-len('.annodbg')] if path.endswith('.annodbg') else path
    head, sep, tail = name.rpartition('.')
    return head if sep and tail and '_' in tail else name


class ManifestError(Exception):
    pass


# parse_annotation_type's annotation extensions, in its order (the first suffix that matches
# decides), with what the loader reads beside each: (extension, the row-diff anchors and fork
# successors beside the GRAPH, the sequence headers). Mirrors metagraph's
# index_annotation_kinds(), which `traverse --index-inventory` prints (annotation_kinds)
ANNOTATION_KINDS = [
    ('.column.annodbg', False, False),
    ('.column_coord.annodbg', False, True),
    ('.brwt_coord.annodbg', False, True),
    ('.row_diff_coord.annodbg', False, True),
    ('.row_diff_brwt_coord.annodbg', False, True),
    ('.row_diff.annodbg', True, False),
    ('.row.annodbg', False, False),
    ('.brwt.annodbg', False, False),
    ('.row_diff_brwt.annodbg', False, False),
    ('.bin_rel_wt.annodbg', False, False),
    ('.flat.annodbg', False, False),
    ('.row_sparse.annodbg', False, False),
    ('.row_diff_flat.annodbg', False, False),
    ('.row_diff_sparse.annodbg', False, False),
    ('.row_diff_disk.annodbg', False, False),
    ('.row_diff_int_disk.annodbg', False, False),
    ('.row_diff_disk_coord.annodbg', False, True),
    ('.rbfish.annodbg', False, False),
    ('.rb_brwt.annodbg', False, False),
    ('.int_brwt.annodbg', False, False),
    ('.row_diff_int_brwt.annodbg', False, False),
]


def _opens(path):
    try:
        with open(path, 'rb'):
            return True
    except OSError:
        return False


# The rule of the graph's derived data, as metagraph states it (kIndexDerivedDataRule;
# `traverse --index-inventory` prints it as derived_rule, and an integration test compares)
DERIVED_RULE = ('the dummy-edge mask (.edgemask) and the Bloom filter (.bloom) are derived data of '
                'the graph, not part of index_fp: they decide which counts are exact and how fast '
                'k-mers are looked up, never what an exact answer is, so adding, removing or '
                'rebuilding one leaves index_fp unchanged, and a manifest lists neither')


def derived_role(path):
    """'graph_mask' or 'graph_bloom' for a file whose base name ends in .edgemask or .bloom
    (the derived data of a graph, whichever graph it belongs to), None for any other; as
    metagraph's index_derived_role. A manifest lists none."""
    name = os.path.basename(path)
    if name.endswith('.edgemask'):
        return 'graph_mask'
    if name.endswith('.bloom'):
        return 'graph_bloom'
    return None


def derived_files(graph):
    """[(path, role, exists, loaded)]: the derived data DBGSuccinct::load reads beside the graph
    listed as |graph| -- <graph without .dbg>.edgemask when it opens, and only then <graph
    without .dbg>.bloom when it exists -- both candidates of a .dbg graph whether they exist or
    not, none for another graph type; as metagraph's index_derived_files. Not part of index_fp
    (DERIVED_RULE)."""
    if not graph.endswith('.dbg'):
        return []
    prefix = graph[:-len('.dbg')]
    mask, bloom = prefix + '.edgemask', prefix + '.bloom'
    mask_loaded = _opens(mask)
    bloom_exists = os.path.exists(bloom)
    return [(mask, 'graph_mask', os.path.exists(mask), mask_loaded),
            (bloom, 'graph_bloom', bloom_exists, mask_loaded and bloom_exists)]


def load_inventory(graph, annotation, coord_mapping=True):
    """[(path, role, required)]: every identity file the server's loaders open for the pair
    listed as |graph| and |annotation|, as metagraph's index_load_inventory derives them (a
    required file is listed whether it exists or not, an optional one only when the loader
    would read it). The graph's derived data (derived_files) is not among them."""
    files = [(graph, 'graph', True)]
    files.append((annotation, 'annotation', True))
    for extension, anchors, coordinates in ANNOTATION_KINDS:
        if not annotation.endswith(extension):
            continue
        if anchors:
            files.append((graph + '.anchors', 'row_diff_anchors', True))
            files.append((graph + '.rd_succ', 'row_diff_fork_succ', True))
        if coordinates and coord_mapping:
            seqs = annotation[:-len(extension)] + '.seqs'
            if os.path.exists(seqs):
                files.append((seqs, 'coord_to_header', False))
        break
    return files


def unloaded_optional(graph, annotation, coord_mapping=True):
    """The optional identity files the inventory derives for the pair -- a coordinate
    annotation's .seqs -- that it does not load (missing beside the listed spelling, or
    coord_mapping off), as metagraph's index_unloaded_optional_files: the server refuses a
    manifest that lists one (by base name), so that index_fp identifies the files it loads.
    (The graph's mask and Bloom filter are refused in any manifest: derived_role.)"""
    candidates = []
    seqs = coordinate_headers(annotation)
    if seqs:
        candidates.append(seqs)
    loaded = [path for path, _, _ in load_inventory(graph, annotation, coord_mapping)]
    return [c for c in candidates if c not in loaded]


def duplicate_base_names(paths):
    """The base names more than one of |paths| has (sorted): the server matches loaded files to
    manifest entries by base name, so it refuses a manifest with any."""
    names = [os.path.basename(p) for p in paths]
    return sorted({n for n in names if names.count(n) > 1})


def sidecars(graph, annotation, coord_mapping=True):
    """The files of the inventory besides the graph and the annotation (the identity files the
    server loads beside them; never the graph's derived mask or Bloom filter)."""
    return [path for path, role, _ in load_inventory(graph, annotation, coord_mapping)
            if role not in ('graph', 'annotation')]


def coordinate_headers(annotation):
    """<base>.<type>.annodbg -> <base>.seqs, where the loader reads the sequence headers of a
    coordinate annotation type (None for another type)."""
    for extension, _, coordinates in ANNOTATION_KINDS:
        if annotation.endswith(extension):
            return annotation[:-len(extension)] + '.seqs' if coordinates else None
    return None


def is_graph_or_annotation(path):
    name = os.path.basename(path)
    return name.endswith('.annodbg') or name.endswith('dbg')


def bundle_files(graph, annotation, extra, *, fail=None, coord_mapping=True):
    files = [path for path, _, _ in load_inventory(graph, annotation, coord_mapping)]
    others = [f for f in extra if is_graph_or_annotation(f)]
    if others:
        (fail or sys.exit)('--extra %s: a manifest describes one graph with one annotation (the '
                           'server refuses one that lists another)' % ', '.join(others))
    derived = [f for f in extra if derived_role(f)]
    if derived:
        (fail or sys.exit)('--extra %s: %s (the server refuses a manifest that lists one)'
                           % (', '.join(derived), DERIVED_RULE))
    unloaded = {os.path.basename(p) for p in unloaded_optional(graph, annotation, coord_mapping)}
    not_loaded = [f for f in extra if os.path.basename(f) in unloaded]
    if not_loaded:
        (fail or sys.exit)('--extra %s: an optional file of the loader inventory that the server '
                           'does not load for this pair (the server refuses a manifest that lists '
                           'it)' % ', '.join(not_loaded))
    files += list(extra)
    missing = [f for f in files if not os.path.isfile(f)]
    if missing:
        (fail or sys.exit)('missing bundle file(s): ' + ', '.join(missing))
    # one entry per file, whatever way it was named
    seen, out = set(), []
    for f in files:
        real = os.path.realpath(f)
        if real not in seen:
            seen.add(real)
            out.append(f)
    return out


def canonical_lines(entries):
    entries = sorted(entries, key=lambda e: e['path'].encode())
    return ''.join(f"{e['path']}\t{e['size']}\t{e['sha256']}\n" for e in entries)


def write(args):
    out = args.output or annotation_base(args.annotation) + '.manifest.json'
    files = bundle_files(args.graph, args.annotation, args.extra,
                         coord_mapping=not args.no_coord_mapping)
    total = sum(os.path.getsize(f) for f in files)
    print(f'{len(files)} files, {total / 1e9:.2f} GB to hash', file=sys.stderr)
    entries, hashed, t0 = [], 0, time.time()
    for f in files:
        size = os.path.getsize(f)

        def progress(done, base=hashed, size=size):
            el = time.time() - t0
            rate = (base + done) / el / 1e6 if el > 0 else 0
            print(f'\r  {os.path.basename(f)}: {done / max(size, 1) * 100:5.1f}%  ({rate:.0f} MB/s)',
                  end='', file=sys.stderr)

        digest = sha256_file(f, progress)
        hashed += size
        print(file=sys.stderr)
        entries.append({'path': os.path.basename(f), 'size': size, 'sha256': digest})
    names = [e['path'] for e in entries]
    if len(set(names)) != len(names):
        sys.exit('two bundle files share a base name: ' + ', '.join(sorted(n for n in names if names.count(n) > 1)))
    entries.sort(key=lambda e: e['path'].encode())
    manifest = {
        'format': 'metagraph-index-manifest',
        'version': 1,
        'index_fp': hashlib.sha256(canonical_lines(entries).encode()).hexdigest(),
        'files': entries,
    }
    if args.name:
        manifest['index_ns'] = args.name
    with open(out, 'w') as f:
        json.dump(manifest, f, indent=1)
        f.write('\n')
    print(f'{out}: index_fp {manifest["index_fp"]} ({time.time() - t0:.0f} s)')


# ---------------------------------------------------------------- batch mode (--server-csv)

def parse_server_csv(path):
    """The lines of a server's graph list as server_query reads them: [(line_no, columns)],
    columns = [name, graph, annotation, manifest, index_ns] ('' for an absent optional one).
    Raises ManifestError for a line the server would refuse."""
    lines = []
    with open(path, encoding='utf-8', newline='') as f:
        text = f.read()
    for no, line in enumerate(text.split('\n'), 1):
        if line == '':
            continue
        cols = line.split(',')
        if len(cols) < 3 or len(cols) > 5:
            raise ManifestError('line %d of %s: %d columns; expected name,graph_path,'
                                'annotation_path[,manifest_path[,index_ns]]'
                                % (no, path, len(cols)))
        cols += [''] * (5 - len(cols))
        if cols[4] and not re.fullmatch(r'[A-Za-z0-9._-]+', cols[4]):
            raise ManifestError('line %d of %s: index_ns %r does not match [A-Za-z0-9._-]+'
                                % (no, path, cols[4]))
        lines.append((no, cols))
    return lines


_DIGEST_LINE = re.compile(r'^([0-9a-fA-F]{64}) ( |\*)(.+)$')


def read_digests(paths):
    """[(digest file, line number, path as written, sha256)] from sha256sum-style files ('#'
    comments, blank lines); one path written twice with different digests is an error."""
    out = []
    for p in paths:
        seen = {}
        with open(p, encoding='utf-8') as f:
            for no, line in enumerate(f, 1):
                line = line.rstrip('\r\n')
                if not line.strip() or line.lstrip().startswith('#'):
                    continue
                m = _DIGEST_LINE.match(line)
                if not m:
                    raise ManifestError('%s line %d: expected "<sha256>  <path>"' % (p, no))
                digest, name = m.group(1).lower(), m.group(3)
                if seen.get(name, digest) != digest:
                    raise ManifestError('%s line %d: %s has two different digests' % (p, no, name))
                seen[name] = digest
                out.append((p, no, name, digest))
    return out


class DigestIndex:
    """Precomputed digests matched to bundle files: by real path (a relative path in a digest
    file is taken from that file's directory), and with |base_names| by a base name that only
    one digest line carries among the lines whose paths name no file here (see the module's
    doc)."""

    def __init__(self, lines, base_names=False):
        self.base_names = base_names
        self.real = {}          # real path -> [(digest, where)]
        self.base = {}          # base name -> [(digest, where)]
        for source, no, name, digest in lines:
            where = '%s line %d (%s)' % (source, no, name)
            path = name if os.path.isabs(name) else os.path.join(
                os.path.dirname(os.path.abspath(source)), name)
            self.real.setdefault(os.path.realpath(path), []).append((digest, where))
            # a line whose path names a file of this machine is that file's digest only: it
            # is never matched to another file by base name, whatever the order of lookups
            if not os.path.exists(path):
                self.base.setdefault(os.path.basename(name), []).append((digest, where))
        self.used = {}          # digest line -> the real path it was matched to
        self.base_matches = []  # (bundle file, digest line), for the report

    def find(self, path):
        """The digest of |path|, None when no digest line names it. Two different digests for
        one file, or one digest line matched by base name to two files, are errors."""
        real = os.path.realpath(path)
        named = self.real.get(real)
        if named:
            if len({d for d, _ in named}) > 1:
                raise ManifestError('two different digests name %s: %s'
                                    % (path, '; '.join(w for _, w in named)))
            for _, where in named:
                self.used[where] = real
            return named[0][0]
        if not self.base_names:
            return None
        cands = self.base.get(os.path.basename(path), [])
        if not cands or len({d for d, _ in cands}) > 1:
            return None         # none, or ambiguous: two different digests
        # a line that names one file (by path, or by base name before) is not another's
        digest, where = cands[0]
        other = self.used.setdefault(where, real)
        if other != real:
            raise ManifestError('the digest of %s would be matched by base name to two '
                                'different files, %s and %s: give the digests by path'
                                % (where, other, real))
        self.base_matches.append((path, where))
        return digest


def batch(args, hasher=sha256_file):
    lines = parse_server_csv(args.server_csv)
    if not lines:
        raise ManifestError('%s lists no graphs' % args.server_csv)
    digests = (DigestIndex(read_digests(args.digests),
                           base_names=getattr(args, 'match_base_names', False))
               if args.digests else None)
    # The unique pairs in order of first appearance, with the first line naming each and the
    # manifest path and index_ns any of their lines give: the server reads one identity per
    # pair (an empty column takes the other line's value, two different ones refuse to start),
    # so a manifest_path given on any line of a pair is where its manifest must go
    pairs = {}
    for no, cols in lines:
        entry = pairs.setdefault((cols[1], cols[2]), {'no': no, 'cols': cols, 'manifest': '',
                                                      'manifest_no': 0, 'ns': '', 'ns_no': 0})
        if cols[3]:
            if (entry['manifest'] and os.path.realpath(entry['manifest'])
                    != os.path.realpath(cols[3])):
                raise ManifestError('lines %d and %d give different manifest paths for one '
                                    '(graph, annotation) pair: %s, %s'
                                    % (entry['manifest_no'], no, entry['manifest'], cols[3]))
            if not entry['manifest']:
                entry['manifest'], entry['manifest_no'] = cols[3], no
        if cols[4]:
            if entry['ns'] and entry['ns'] != cols[4]:
                raise ManifestError('lines %d and %d give different index_ns for one '
                                    '(graph, annotation) pair: %s, %s'
                                    % (entry['ns_no'], no, entry['ns'], cols[4]))
            entry['ns'], entry['ns_no'] = cols[4], no
    plan = []          # (pair, line_no, index_ns, files, out)
    outputs = {}
    for (graph, annotation), entry in pairs.items():
        no, cols = entry['no'], entry['cols']

        def fail(msg, no=no):
            raise ManifestError('line %d: %s' % (no, msg))
        files = bundle_files(graph, annotation, [], fail=fail,
                             coord_mapping=not getattr(args, 'no_coord_mapping', False))
        if entry['manifest']:
            out = entry['manifest']
        elif args.out_dir:
            out = os.path.join(args.out_dir, '%s.%s.manifest.json'
                               % (cols[0], os.path.basename(annotation_base(annotation))))
        else:
            out = annotation_base(annotation) + '.manifest.json'
        key = os.path.realpath(out)
        if key in outputs:
            raise ManifestError('lines %d and %d would both write %s' % (outputs[key], no, out))
        outputs[key] = no
        if os.path.exists(out) and not args.force:
            raise ManifestError('line %d: %s exists (--force replaces it)' % (no, out))
        names = [os.path.basename(f) for f in files]
        if len(set(names)) != len(names):
            raise ManifestError('line %d: two bundle files share a base name' % no)
        plan.append(((graph, annotation), no, entry['ns'], files, out))
    # every distinct file once (by real path): a precomputed digest, else hashed
    distinct = {}
    for _, _, _, files, _ in plan:
        for f in files:
            distinct.setdefault(os.path.realpath(f), f)
    known, todo = {}, []
    for real, f in distinct.items():
        matches = len(digests.base_matches) if digests else 0
        d = digests.find(f) if digests else None
        if d is not None:
            by_base = len(digests.base_matches) > matches
            known[real] = (d, 'precomputed-by-base-name' if by_base else 'precomputed')
        elif args.digests_only:
            raise ManifestError('no precomputed digest for %s (--digests-only)' % f)
        else:
            todo.append((real, f))
    total = sum(os.path.getsize(f) for _, f in todo)
    print('%d pairs, %d distinct files: %d precomputed, %d to hash (%.2f GB, %d streams)'
          % (len(plan), len(distinct), len(known), len(todo), total / 1e9, args.jobs),
          file=sys.stderr)
    for f, where in (digests.base_matches if digests else []):
        print('  by base name: %s <- %s' % (f, where), file=sys.stderr)
    lock = threading.Lock()
    done = [0]
    t0 = time.time()

    def hash_one(item):
        real, f = item
        size = os.path.getsize(f)
        d = hasher(f)
        with lock:
            done[0] += size
            el = time.time() - t0
            print('  %s  %s  (%.0f MB/s overall)' % (d[:12], f, done[0] / max(el, 1e-9) / 1e6),
                  file=sys.stderr)
        return real, d

    with ThreadPoolExecutor(max(1, args.jobs)) as ex:
        for real, d in ex.map(hash_one, todo):
            known[real] = (d, 'hashed')
    elapsed = time.time() - t0
    written = {}
    for pair, no, ns, files, out in plan:
        entries = []
        for f in files:
            d, how = known[os.path.realpath(f)]
            entries.append({'path': os.path.basename(f), 'size': os.path.getsize(f),
                            'sha256': d, 'digest': how})
        entries.sort(key=lambda e: e['path'].encode())
        manifest = {
            'format': 'metagraph-index-manifest',
            'version': 1,
            'index_fp': hashlib.sha256(canonical_lines(entries).encode()).hexdigest(),
            'files': entries,
        }
        if ns:
            manifest['index_ns'] = ns
        d = os.path.dirname(os.path.abspath(out))
        os.makedirs(d, exist_ok=True)
        tmp = out + '.tmp'
        with open(tmp, 'w') as f:
            json.dump(manifest, f, indent=1)
            f.write('\n')
        os.replace(tmp, out)
        written[pair] = out
        print('%s: (%s, %s) index_fp %s' % (out, pair[0], pair[1], manifest['index_fp']))
    if args.write_csv:
        with open(args.write_csv, 'w', encoding='utf-8', newline='') as f:
            for no, cols in lines:
                # every line of a pair names the manifest written for it (a line that named
                # it keeps its own spelling)
                row = cols[:3] + [cols[3] or written[(cols[1], cols[2])]]
                if cols[4]:
                    row.append(cols[4])
                f.write(','.join(row) + '\n')
        print('%s: the graph list with its manifests' % args.write_csv)
    print('%.2f GB hashed in %.0f s (%.0f MB/s)' % (total / 1e9, elapsed,
                                                    total / max(elapsed, 1e-9) / 1e6),
          file=sys.stderr)
    return written


def verify(path, root=None):
    with open(path) as f:
        manifest = json.load(f)
    root = root or os.path.dirname(os.path.abspath(path))
    problems = []
    # the server matches loaded files by base name: entries sharing one (a manifest written for
    # a directory of bundles) would lend one index_fp to every bundle, and it refuses them
    for name in duplicate_base_names([e['path'] for e in manifest['files']]):
        problems.append(f'base name listed more than once: {name} (a manifest describes one '
                        f'bundle; the server refuses it)')
    # the graph's derived data is not part of index_fp: the server refuses a manifest that
    # lists it, however its digests compare
    for e in manifest['files']:
        if derived_role(e['path']):
            problems.append(f'derived data listed: {e["path"]} ({DERIVED_RULE}; the server '
                            f'refuses this manifest: write it again without the entry, which '
                            f'changes its index_fp)')
    for e in manifest['files']:
        p = os.path.join(root, e['path'])
        if not os.path.isfile(p):
            problems.append(f'missing: {e["path"]}')
            continue
        if os.path.getsize(p) != e['size']:
            problems.append(f'size differs: {e["path"]}')
            continue
        if sha256_file(p) != e['sha256']:
            problems.append(f'content differs: {e["path"]}')
    fp = hashlib.sha256(canonical_lines(manifest['files']).encode()).hexdigest()
    if manifest.get('index_fp') != fp:
        problems.append(f'index_fp {manifest.get("index_fp")} is not the digest of the file list ({fp})')
    if problems:
        print('\n'.join(problems))
        sys.exit(1)
    print(f'{path}: OK, index_fp {fp}')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('-i', '--graph')
    ap.add_argument('-a', '--annotation')
    ap.add_argument('-o', '--output', help='default: <annotation base>.manifest.json')
    ap.add_argument('--name', help='index_ns to record (the server takes it from --index-name)')
    ap.add_argument('--extra', action='append', default=[], help='a further bundle file')
    ap.add_argument('--verify', metavar='MANIFEST')
    ap.add_argument('--root', help='--verify: the directory holding the bundle (default: next to the manifest)')
    ap.add_argument('--server-csv', metavar='GRAPHS.csv',
                    help='batch mode: one manifest per (graph, annotation) pair of a server list')
    ap.add_argument('--jobs', type=int, default=4, help='batch: parallel hashing streams [4]')
    ap.add_argument('--digests', action='append', default=[], metavar='FILE',
                    help='batch: precomputed "<sha256>  <path>" lines (repeatable)')
    ap.add_argument('--digests-only', action='store_true',
                    help='batch: hash nothing; a file without a precomputed digest is an error')
    ap.add_argument('--match-base-names', action='store_true',
                    help='batch: also match digests by a base name only one digest line has '
                         '(digest files whose paths are not this machine\'s); each match is listed')
    ap.add_argument('--out-dir', help='batch: where manifests without a CSV column go')
    ap.add_argument('--write-csv', metavar='OUT.csv', help='batch: the list with its manifests')
    ap.add_argument('--force', action='store_true', help='batch: replace existing manifests')
    ap.add_argument('--no-coord-mapping', action='store_true',
                    help='the server runs with --no-coord-mapping: it does not load the .seqs')
    ap.add_argument('--inventory', action='store_true',
                    help='with -i and -a: print the loader dependency inventory (JSON: the '
                         'identity files, and the derived mask and Bloom filter no manifest '
                         'lists) and exit')
    args = ap.parse_args()
    if args.verify:
        verify(args.verify, args.root)
    elif args.server_csv:
        try:
            batch(args)
        except ManifestError as e:
            sys.exit('index_manifest.py: %s' % e)
    elif args.graph and args.annotation and args.inventory:
        # as `metagraph traverse --index-inventory`: the identity files a manifest lists, and the
        # graph's derived data the loader reads beside them, which no manifest lists
        print(json.dumps({
            'files': [{'path': p, 'role': r, 'required': q, 'exists': os.path.exists(p)}
                      for p, r, q in load_inventory(args.graph, args.annotation,
                                                    not args.no_coord_mapping)],
            'derived': [{'path': p, 'role': r, 'exists': e, 'loaded': l}
                        for p, r, e, l in derived_files(args.graph)],
            'derived_rule': DERIVED_RULE}, indent=1))
    elif args.graph and args.annotation:
        write(args)
    else:
        ap.error('give -i and -a, --server-csv, or --verify')


if __name__ == '__main__':
    main()
