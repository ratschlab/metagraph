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
    index_manifest.py --verify MANIFEST          # re-hash every listed file and compare

The bundle is the graph and the annotation plus the sidecars found next to them with a known
suffix (graph: .anchors, .rd_succ, .edgemask, .weights; annotation base: .seqs, .coords).
Anything else that belongs to the bundle is added with --extra. Paths in the manifest are
the files' BASE NAMES, so the identity does not depend on where the bundle or the manifest
lives (the server matches loaded files by base name too); the base names in one bundle must be
distinct. --verify looks for the files next to the manifest, or under --root. Stdlib only.
"""
import argparse
import hashlib
import json
import os
import sys
import time

GRAPH_SIDECARS = ('.anchors', '.rd_succ', '.edgemask', '.weights')
ANNOTATION_SIDECARS = ('.seqs', '.coords')


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


def bundle_files(graph, annotation, extra):
    files = [graph, annotation]
    files += [graph + s for s in GRAPH_SIDECARS if os.path.isfile(graph + s)]
    base = annotation_base(annotation)
    files += [base + s for s in ANNOTATION_SIDECARS if os.path.isfile(base + s)]
    files += list(extra)
    missing = [f for f in files if not os.path.isfile(f)]
    if missing:
        sys.exit('missing bundle file(s): ' + ', '.join(missing))
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
    files = bundle_files(args.graph, args.annotation, args.extra)
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


def verify(path, root=None):
    with open(path) as f:
        manifest = json.load(f)
    root = root or os.path.dirname(os.path.abspath(path))
    problems = []
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
    args = ap.parse_args()
    if args.verify:
        verify(args.verify, args.root)
    elif args.graph and args.annotation:
        write(args)
    else:
        ap.error('give -i and -a, or --verify')


if __name__ == '__main__':
    main()
