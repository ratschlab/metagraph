"""pandas views of a graphlet. pandas is imported here, lazily, and nowhere else in the
package: `import metagraph.traverse` must work without it.

frames(g) -> {name: DataFrame} with one row per labels / segments / runs / walks /
splits / claims / events (and presence runs in annotate mode). Label columns hold refs
(the stable join key) with the name beside them; ids are the format's ordinals.
"""

from . import derive
from .model import ARM_SIDES

__all__ = ['frames', 'records']


def records(g):
    """The rows frames() turns into DataFrames, as plain lists of dicts (no pandas)."""
    out = {k: [] for k in ('labels', 'segments', 'runs', 'walks', 'splits', 'claims',
                           'events', 'presence')}
    for l in g.labels:
        out['labels'].append({'id': l.id, 'ref': l.ref, 'name': l.name, 'kind': l.kind,
                              'column': l.column, 'seq_id': l.seq_id})
    for side in ARM_SIDES:
        a = g.arms.get(side)
        if a is None:
            continue
        for s in a.segments:
            out['segments'].append({
                'arm': side, 'segment': s.id, 'parents': list(s.parents),
                'from_bp': s.from_bp, 'to_bp': s.end_bp, 'length_bp': s.length_bp,
                'labels_in': len(s.entry), 'labels_in_total': s.entry_total,
                'labels_out': len(s.end), 'leaf': s.leaf is not None,
                'merge': len(s.parents) > 1, 'split': s.split})
            for p in s.presence:
                out['presence'].append({'arm': side, 'segment': s.id, 'from_bp': p.from_bp,
                                        'to_bp': p.to_bp, 'labels': len(p.labels),
                                        'labels_total': p.total})
            for ev in s.events:
                out['events'].append({'arm': side, 'segment': s.id, 'at_bp': ev.at_bp,
                                      'type': ev.type, 'char': ev.char,
                                      'from': ev.from_label, 'to': ev.to_label,
                                      'cost': ev.cost})
        for r in a.runs:
            out['runs'].append({
                'arm': side, 'run': r.id, 'segment': r.segment, 'ref': g.labels[r.label].ref,
                'name': g.labels[r.label].name, 'from_bp': r.from_bp, 'to_bp': r.to_bp,
                'end': r.end, 'reason': r.reason, 'qualifier': r.qualifier,
                'entered_by': r.entered_by, 'route_bp': r.route_bp,
                'evidence_from': derive.evidence(a, r)[1], 'loss': r.loss,
                'branches': r.branches})
        for w in g.walks(a, by='id'):
            out['walks'].append({
                'arm': side, 'walk': w.path_id, 'leaf': w.leaf, 'length_bp': w.length_bp,
                'path_reason': w.path_reason, 'n_alive': w.n_alive,
                'n_full': len(w.labels_full), 'complete': w.complete,
                'beyond_certified_bp': w.beyond_certified_bp})
        for sp in derive.splits(a, g.mode):
            out['splits'].append({'arm': side, 'at_bp': sp.at_bp, 'segment': sp.segment,
                                  'children': list(sp.children),
                                  'kind': 'ambiguous' if sp.ambiguous else 'divergence',
                                  'labels_before': sp.labels_before})
        for c in g.claims(a, strict=False):
            d = c.as_dict()
            d['ref'], d['name'] = c.label.ref, c.label.name
            del d['label']
            d['from_label'] = c.from_label.ref if c.from_label else None
            out['claims'].append(d)
    return out


def frames(g):
    import pandas as pd   # lazy: the package itself never needs pandas
    return {k: pd.DataFrame(v) for k, v in records(g).items()}

