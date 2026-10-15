"""The oracle run's report: oracle_report.json (everything) and oracle_report.md (for people).

The markdown states the scope (what phase 1 covers and what it does not), the outcome per
config, every disagreement grouped by check (the request, the oracle's set and the engine's
answer), the exploratory differences, and the corruptions the comparison missed.
"""

import json

COVERED = [
    '/pattern, patterns of at most k bases: dna, iupac, protein (genetic codes 1, 2, 4, 6, 11, '
    '27, 31), strands both/forward/reverse, scopes any_offset and suffix, max_mismatches 0-3 with '
    'by_mismatches levels and substitutions, modes count/all_or_count/partial, labels "all" with '
    'record placement (BASIC) and labels of both strands (CANONICAL, PRIMARY)',
    'graph modes BASIC (k 10, 11, 15), CANONICAL and wrapped PRIMARY (k 10 and 11: even k with '
    'k-mers equal to their reverse complement); masked graphs, and BASIC k 11 without its mask at '
    'max_checked_entries 50 and 0, every bounds count held to the SPEC\'s U (the masked count plus '
    'the source dummies holding the pattern, per level)',
    '/resolve with explicit column and header labels, support kmer and trace, on BASIC',
    'budget probes (max_steps, stop_at_threshold): invariants of legitimately incomplete answers',
    '/traverse (phase 2): annotate\'s exhaustive trie and constrain\'s walks under keep as complete '
    'sets (walks, leaves, label ends with reason, branches and loss, node label sets, blocked and '
    'hairpin events, seed and dropped labels); forbid and constant costs, switch_on loss/any, the '
    'branch allowance, the quorum, min_live_labels, the split limit, trace support, hairpins skip '
    'and follow, edge_reuse(_rc), rejoined_seed; cap probes (max_steps, max_paths, '
    'max_live_paths) judged up to complete_to_bp; merge: an arm that merged nothing held to the '
    'keep walker, an arm that merged held to it before its first merge and judged after it on the '
    'displays through first parents, the merge structure (parents, labels_via_parent, the '
    'majority parent, reconverge events), each leaf label\'s own route and displayed support; '
    'tree and summary against full',
]
NOT_COVERED = [
    '/traverse: table and hll costs, the graphlet\'s MGT body, output.coordinates, request '
    'budgets and attempts, beam pruning, growth bins; after a merge, the united history\'s label '
    'states (which labels reach a node, with which loss and branches), its run ends and events, '
    'its needed budgets and which walks it continues; ties of the majority parent (arrival '
    'order)', 'the row of a CANONICAL k-mer listed without its reverse complement, when it is '
    'not the canonical one: held only to row < node - 1', 'patterns longer than k (anchors, long_search paths and supported '
    'paths)', 'predicates and motif predicates', 'multi-graph servers and in_ram',
    'DNA5 graphs (refused by this build: alphabet_untested)', 'CANONICAL and PRIMARY graphs '
    'without their mask', 'the lower bound of a bounds count beyond lower <= truth, and the '
    'estimate beyond lower <= estimate <= upper', 'time stops and the deadline',
    'the low-complexity note', 'work counters beyond their stated sums',
    '/resolve discover, select and run_format "rle"', 'the information floor (lowered to 0)',
]


def engine_provenance(path):
    """What was actually run (verification M1: "--version" alone hid a build of another commit):
    the binary's real path, size, sha256 and modification time, its --version, and the git HEAD of
    the tree that holds it with how many of that tree's metagraph/src files are uncommitted."""
    import datetime
    import hashlib
    import os
    import subprocess
    out = {'path': path}
    try:
        real = os.path.realpath(path)
        h = hashlib.sha256()
        with open(real, 'rb') as f:
            for block in iter(lambda: f.read(1 << 20), b''):
                h.update(block)
        st = os.stat(real)
        out.update({'realpath': real, 'size': st.st_size, 'sha256': h.hexdigest(),
                    'mtime_utc': datetime.datetime.fromtimestamp(
                        st.st_mtime, datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')})
    except OSError as e:
        out['error'] = str(e)
        return out
    try:
        out['version'] = subprocess.run([real, '--version'], capture_output=True, text=True,
                                        timeout=60).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        out['version'] = None
    d = os.path.dirname(real)

    def git(*args):
        res = subprocess.run(['git', '-C', d] + list(args), capture_output=True, text=True)
        return res.stdout.strip() if res.returncode == 0 else None
    top = git('rev-parse', '--show-toplevel')
    out['build_tree'] = top
    out['build_head'] = git('rev-parse', 'HEAD') if top else None
    if top:
        # the pathspec from the tree's top: relative to the binary's directory (metagraph/build)
        # it named metagraph/build/metagraph/src and always counted 0
        status = git('status', '--porcelain', '--', ':/metagraph/src')
        out['build_src_uncommitted'] = len(status.splitlines()) if status is not None else None
    return out


def _compact(value, limit=3000):
    text = json.dumps(value, sort_keys=True)
    return text if len(text) <= limit else text[:limit] + ' ...'


def oracle_summary(exp, limit=60):
    if exp is None or getattr(exp, 'pattern', None) is None or exp.pattern.slot:
        return None
    return {'total': exp.counts['total'], 'suffix': exp.counts['suffix'],
            'by_offset': exp.counts['by_offset'], 'by_orient': exp.counts['by_orient'],
            'by_mismatches': exp.counts['by_mismatches'],
            'contexts': [[c.orientation, c.kmer, c.offset, c.level] for c in exp.ctxs[:limit]],
            'contexts_listed': min(limit, len(exp.ctxs))}


def engine_excerpt(entry, limit=60):
    if not isinstance(entry, dict):
        return entry
    out = {k: v for k, v in entry.items() if k not in ('timing',)}
    if isinstance(out.get('results'), list) and len(out['results']) > limit:
        out['results'] = out['results'][:limit] + [f'... {len(entry["results"]) - limit} more']
    return out


def write(path_json, path_md, report):
    with open(path_json, 'w') as f:
        json.dump(report, f, indent=1, sort_keys=True)
        f.write('\n')
    with open(path_md, 'w') as f:
        f.write(markdown(report))


def markdown(report):
    lines = ['# Oracle suite report (lane O: /pattern, /resolve, /traverse)', '']
    eng = report['engine']
    lines.append(f"Engine: `{eng.get('realpath', eng['path'])}`, sha256 `{eng.get('sha256', '?')}`, "
                 f"modified {eng.get('mtime_utc', '?')}, {eng.get('version') or '?'}; build tree "
                 f"`{eng.get('build_tree')}` at HEAD `{eng.get('build_head')}` "
                 f"({eng.get('build_src_uncommitted', '?')} uncommitted files under metagraph/src). "
                 f"Seed {report['seed']}. {report['started_utc']} to {report['ended_utc']}.")
    lines.append('')
    s = report['summary']
    lines.append(f"**Verdict: {s['verdict']}.** {s['entries']} answers compared "
                 f"({s['requests']} requests): {s['pass']} pass, {s['fail']} fail, "
                 f"{s['inconclusive']} inconclusive (budget and cap probes, merge arms after "
                 f"their first merge: completeness not judged, invariants judged), {s['exploratory']} with "
                 f"an exploratory difference. Mutations of the stated mutation set: "
                 f"{report['mutations']['caught']} of {report['mutations']['applied']} caught, "
                 f"{len(report['mutations']['missed'])} missed.")
    lines.append('')
    lines.append('## Scope')
    lines.append('')
    lines.append('Covered:')
    lines += [f'- {c}' for c in COVERED]
    lines.append('')
    lines.append('Not covered (stated, not tested):')
    lines += [f'- {c}' for c in NOT_COVERED]
    lines.append('')
    lines.append('## Per config')
    lines.append('')
    lines.append('| config | graph | k | counting | requests | answers | pass | fail | '
                 'inconclusive | exploratory |')
    lines.append('|---|---|---|---|---|---|---|---|---|---|')
    for name, c in report['configs'].items():
        o = c['outcomes']
        lines.append(f"| {name} | {c.get('graph_mode')} | {c.get('k')} | {c.get('counting')} | "
                     f"{c['requests']} | {c['entries']} | {o.get('pass', 0)} | {o.get('fail', 0)} | "
                     f"{o.get('inconclusive', 0)} | {o.get('exploratory', 0)} |")
    lines.append('')
    lines.append('## Disagreements')
    lines.append('')
    dis = report['disagreements']
    if not dis:
        lines.append('None: every answer equals the oracle (complete sets, exact counts).')
    groups = {}
    for d in dis:
        key = (d['config'], tuple(sorted({c['check'].split('[')[0] for c in d['failed']})))
        groups.setdefault(key, []).append(d)
    for (config, checks), items in sorted(groups.items()):
        lines.append(f"### {config}: {', '.join(checks)} ({len(items)} answers)")
        lines.append('')
        for d in items[:3]:
            lines.append(f"- `{d['request_id']}` pattern {d['position']}: request "
                         f"`{_compact(d['request'], 1500)}`")
            for c in d['failed'][:6]:
                lines.append(f"  - `{c['check']}`: expected `{_compact(c.get('expected'), 800)}`, "
                             f"got `{_compact(c.get('got'), 800)}`")
            if d.get('oracle'):
                lines.append(f"  - oracle: `{_compact(d['oracle'], 1500)}`")
        if len(items) > 3:
            lines.append(f'- ... {len(items) - 3} more in oracle_report.json')
        lines.append('')
    expl = report.get('exploratory', [])
    lines.append('## Exploratory differences')
    lines.append('')
    if not expl:
        lines.append('None.')
    for d in expl[:10]:
        lines.append(f"- `{d['request_id']}`: {d['note']}; expected `{_compact(d.get('expected'), 800)}`, "
                     f"got `{_compact(d.get('got'), 800)}`")
    lines.append('')
    lines.append('## Mutations the comparison missed')
    lines.append('')
    m = report['mutations']
    lines.append('| mutation | applied | caught |')
    lines.append('|---|---|---|')
    for name, v in sorted(m['by_mutation'].items()):
        lines.append(f"| {name} | {v['applied']} | {v['caught']} |")
    lines.append('')
    if not m['missed']:
        lines.append('None missed.')
    for x in m['missed'][:30]:
        lines.append(f"- `{x['request_id']}` pattern {x.get('position')}: `{x['mutation']}` "
                     f"(families {', '.join(x['families'])})")
    lines.append('')
    return '\n'.join(lines) + '\n'
