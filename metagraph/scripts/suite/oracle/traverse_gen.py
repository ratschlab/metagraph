"""Deterministic /traverse requests (ORACLE-DESIGN.md section 10.4).

- random seeds: windows of the records (k to 40 bp, at random and at the features), radius 0-40, each
  direction, annotate and constrain (forbid; explicit labels, some not supporting the seed, or a
  derived set), exhaustive or on_reconverge keep with each branch allowance, trace on BASIC, hairpins
  followed or skipped;
- the traversal fixture's rule cases, each in the setting that permits and the one that forbids;
- cap probes (max_steps, max_paths, max_live_paths), whose answers are legitimately incomplete;
- on_reconverge merge (the default outside the preset), judged on its invariants.
"""

import random

from . import requests_gen

RADII = [1, 3, 8, 15, 25, 40]   # the server refuses 0 (range [1, ...])
RULE_STRATEGIES = {
    'divergence': [{'branching': {'on_reconverge': 'keep'}},
                   {'branching': {'on_reconverge': 'keep', 'max_splits_per_path': 0}},
                   {'exhaustive': True}],
    'ambiguous': [{'branching': {'on_reconverge': 'keep', 'max_label_branches': 0}},
                  {'branching': {'on_reconverge': 'keep', 'max_label_branches': 1}},
                  {'exhaustive': True}],
    'cycle': [{'branching': {'on_reconverge': 'keep', 'max_label_branches': b}} for b in (0, 1, 2)]
             + [{'exhaustive': True}, {'exhaustive': True, 'labels': {'mode': 'annotate'}}],
    'bubble': [{'branching': {'on_reconverge': 'keep'}}, {}, {'exhaustive': True}],
    # sw1 -> sw2 -> sw3: forbid loses sw1; budget 1 takes one switch and ends sw2 with
    # loss_budget (needed 2); budget 2 takes both; switch_on any enters early
    'transition': [{'branching': {'on_reconverge': 'keep'}},
                   {'branching': {'on_reconverge': 'keep'}, 'labels': {
                       'change_cost': {'model': 'constant', 'value': 1}, 'loss_budget': 1}},
                   {'branching': {'on_reconverge': 'keep'}, 'labels': {
                       'change_cost': {'model': 'constant', 'value': 1}, 'loss_budget': 2}},
                   {'branching': {'on_reconverge': 'keep'}, 'labels': {
                       'change_cost': {'model': 'constant', 'value': 1}, 'loss_budget': 2,
                       'switch_on': 'any'}}],
    'fork': [{'branching': {'on_reconverge': 'keep', 'max_label_branches': b}, 'labels': {
                 'change_cost': {'model': 'constant', 'value': 1}, 'loss_budget': 2}}
             for b in (0, 1)]
            + [{'branching': {'on_reconverge': 'keep'}, 'labels': {
                'change_cost': {'model': 'constant', 'value': 1}, 'loss_budget': 2,
                'switch_on': 'any'}},
               {'branching': {'on_reconverge': 'keep'}}],
    'quorum': [{'branching': {'on_reconverge': 'keep'}},
               {'branching': {'on_reconverge': 'keep', 'min_successor_labels': 2}},
               {'branching': {'on_reconverge': 'keep', 'min_successor_fraction': 0.5}},
               {'branching': {'on_reconverge': 'keep'}, 'bounds': {'min_live_labels': 2}}],
    'circle': [{'branching': {'on_reconverge': 'keep'}}, {'exhaustive': True},
               {'exhaustive': True, 'labels': {'mode': 'annotate'}}],
    'homopolymer': [{'exhaustive': True}, {'branching': {'on_reconverge': 'keep',
                                                         'max_label_branches': 1}},
                    {'exhaustive': True, 'labels': {'mode': 'annotate'}}],
    'hairpin': [{'exhaustive': True}, {'exhaustive': True, 'branching': {'hairpins': 'follow'}},
                {'exhaustive': True, 'labels': {'mode': 'annotate'},
                 'branching': {'hairpins': 'follow'}}],
    'trace': [{'exhaustive': True}, {'exhaustive': True, 'support': 'trace'},
              {'branching': {'on_reconverge': 'keep'}, 'support': 'trace'},
              {'branching': {'on_reconverge': 'keep'}}],
    'trace_break': [{'exhaustive': True}, {'exhaustive': True, 'support': 'trace'}],
}
# the radius of a case (the transition's chain of two switches is 100 bases long)
CASE_RADIUS = {'transition': 110}


def _merge(base, extra):
    out = {}
    for k in set(base) | set(extra):
        a, b = base.get(k), extra.get(k)
        if isinstance(a, dict) and isinstance(b, dict):
            out[k] = _merge(a, b)
        else:
            out[k] = b if k in extra else a
    return out


class TraverseGen:
    def __init__(self, index, walker, seed, config_name, cases=None):
        self.index = index
        self.w = walker
        self.k = index.k
        self.rng = random.Random(f'{seed}:{config_name}:traverse')
        self.pg = requests_gen.Gen(index, seed, config_name + ':traverse')
        self.pg.rng = self.rng
        self.cases = cases
        self.config = config_name
        self.seed = seed
        self.n = 0
        self.header = index.mode == 'basic' and bool(index.seq_ids)

    def next_id(self):
        self.n += 1
        return f'T-{self.config}-{self.n:04d}'

    def label_names(self, names):
        """A case's label names as this index serves them: headers on BASIC with record mapping,
        else the records' columns."""
        if self.header:
            return names
        out = []
        for n in names:
            col = n if n in self.index.columns else next(
                r.column for r in self.index.records if r.name == n)
            if col not in out:
                out.append(col)
        return out

    def _labels_for(self, source, seed):
        rng = self.rng
        kind = 'header' if self.header and rng.random() < 0.7 else 'column'
        pool = [l.name for l in self.w.labels_of_kind(kind)]
        if source is None:
            return None
        own = source.name if kind == 'header' else source.column
        names = [own] + rng.sample(pool, rng.randint(0, 3))
        names = list(dict.fromkeys(names))
        rng.shuffle(names)
        return names

    def random_requests(self, count):
        rng = self.rng
        out = []
        for _ in range(count):
            strategy = {'output': {'detail': 'full'}}
            mode = 'annotate' if rng.random() < 0.3 else 'constrain'
            radius = rng.choice(RADII)
            strategy['bounds'] = {'max_extension_bp': radius}
            d = rng.choice(['both', 'both', 'left', 'right'])
            if d != 'both':
                strategy['direction'] = d
            roll = rng.random()
            if mode == 'annotate' or roll < 0.45:
                strategy['exhaustive'] = True
            elif roll < 0.85:
                strategy['branching'] = {'on_reconverge': 'keep',
                                         'max_label_branches': rng.choice([0, 1, 2, 'unlimited'])}
            # else: the defaults (merge), judged on invariants
            if self.index.mode != 'basic' and rng.random() < 0.3:
                strategy.setdefault('branching', {})['hairpins'] = 'follow'
            if mode == 'annotate':
                strategy['labels'] = {'mode': 'annotate'}
            seeds = []
            for _ in range(rng.choice([1, 1, 2])):
                window, source = self.pg.sourced_window(rng.randint(self.k, 40), None)
                seed = {'sequence': window}
                if mode == 'constrain' and rng.random() < 0.7 and source is not None:
                    seed['labels'] = self._labels_for(source, window)
                seeds.append(seed)
            # trace with merging is refused (400): "a coordinate trace cannot justify" one
            # parent's evidence per label
            merging = not strategy.get('exhaustive') and (strategy.get('branching') or {}).get(
                'on_reconverge') != 'keep'
            if mode == 'constrain' and rng.random() < 0.5 and all('labels' in s for s in seeds):
                self._constant_cost(strategy, seeds)
            if (mode == 'constrain' and self.header and not merging and rng.random() < 0.2
                    and all('labels' not in s or all(n in [r.name for r in self.index.indexed]
                                                    for n in s['labels']) for s in seeds)):
                strategy['support'] = 'trace'
            out.append({'id': self.next_id(), 'route': 'traverse', 'class': 'ample',
                        'request': {'seeds': seeds, 'strategy': strategy}})
        return out

    def _constant_cost(self, strategy, seeds):
        """A constant change cost drawn per request (verification M3): its value, budget and
        switch_on, and up to three extra labels of the seed labels' kind that no seed names;
        an extra above the budget is unreachable, a refusal the oracle predicts (SPEC 5)."""
        rng = self.rng
        value = rng.choice([0.5, 1, 1, 2])
        budget = rng.choice([0, 0.5, 1, 1, 2, 3]) if rng.random() < 0.15 else \
            rng.choice([value, value, 2 * value, 3 * value])
        lab = strategy.setdefault('labels', {})
        lab['change_cost'] = {'model': 'constant', 'value': value}
        lab['loss_budget'] = budget
        if rng.random() < 0.4:
            lab['switch_on'] = 'any'
        named = {n for s in seeds for n in s.get('labels') or []}
        kinds = {self.w.label(n).kind for n in named if self.w.label(n) is not None}
        kind = 'header' if kinds == {'header'} else 'column'
        pool = [l.name for l in self.w.labels_of_kind(kind) if l.name not in named]
        if pool and rng.random() < 0.8:
            lab['extra'] = rng.sample(pool, min(len(pool), rng.randint(1, 3)))

    def cap_requests(self, count):
        """Cap probes: the incomplete member of the triplets (SPEC 6.10)."""
        rng = self.rng
        out = []
        for _ in range(count):
            mode = 'annotate' if rng.random() < 0.4 else 'constrain'
            cap = rng.choice(['max_steps', 'max_paths', 'max_live_paths'])
            value = {'max_steps': rng.randint(2, 40), 'max_paths': rng.randint(1, 3),
                     'max_live_paths': rng.randint(1, 2)}[cap]
            strategy = {'exhaustive': True, 'output': {'detail': 'full'},
                        'bounds': {'max_extension_bp': rng.choice([15, 25, 40]), cap: value}}
            if mode == 'annotate':
                strategy['labels'] = {'mode': 'annotate'}
            window, source = self.pg.sourced_window(rng.randint(self.k, 30), True)
            seed = {'sequence': window}
            if mode == 'constrain' and source is not None and rng.random() < 0.5:
                seed['labels'] = self._labels_for(source, window)
            out.append({'id': self.next_id(), 'route': 'traverse', 'class': 'budget',
                        'request': {'seeds': [seed], 'strategy': strategy}})
        return out

    def rule_requests(self):
        """Every rule case of the traversal fixture in its permitting and forbidding settings."""
        out = []
        if not self.cases:
            return out
        for case, (seed, labels, extra) in self.cases.items():
            if case.startswith('trace') and not self.header:
                continue
            for strat in RULE_STRATEGIES[case]:
                strategy = _merge({'output': {'detail': 'full'},
                                   'bounds': {'max_extension_bp': CASE_RADIUS.get(case, 45)}},
                                  strat)
                annotate = (strategy.get('labels') or {}).get('mode') == 'annotate'
                req_seed = {'sequence': seed}
                if not annotate:
                    req_seed['labels'] = self.label_names(labels)
                    model = ((strategy.get('labels') or {}).get('change_cost') or {}).get('model')
                    # an extra label is a switch target: unreachable (a 400) under forbid
                    if extra and model not in (None, 'forbid'):
                        strategy = _merge(strategy, {'labels': {
                            'extra': [n for n in self.label_names(extra)
                                      if n not in req_seed['labels']]}})
                    if strategy.get('support') == 'trace' and not self.header:
                        continue
                out.append({'id': self.next_id(), 'route': 'traverse', 'class': 'rule',
                            'case': case, 'request': {'seeds': [req_seed],
                                                      'strategy': strategy}})
        return out


    def chain_requests(self, count):
        """Constant-cost requests on the transition chain sw1 -> sw2 -> sw3 at random seeds,
        values, budgets, switch_on and branching (verification round 2). The random requests
        rarely meet a lineage that has paid one switch and needs another; here a seed in sw1
        walking right, or in sw3 walking left, has two switches ahead, so a budget below twice the
        value ends the second label loss_budget (needed_budget twice the value) and one of at
        least twice the value takes both. Its own random stream: the other requests stay as
        they were."""
        if not self.cases or 'transition' not in self.cases:
            return []
        rng = random.Random(f'{self.seed}:{self.config}:chains')
        recs = {r.name: r for r in self.index.records}
        chain = ['sw1', 'sw2', 'sw3']
        out = []
        for _ in range(count):
            start = rng.choice(['sw1', 'sw3'])
            arm = 'right' if start == 'sw1' else 'left'
            r = recs[start]
            L = rng.randint(self.k, 25)
            p = rng.randrange(0, len(r.seq) - L + 1)
            value = rng.choice([0.5, 1, 2, 3])
            labels = {'change_cost': {'model': 'constant', 'value': value},
                      'loss_budget': value * rng.choice([1, 1, 1.5, 1.5, 2, 3])}
            if rng.random() < 0.4:
                labels['switch_on'] = 'any'
            names = self.label_names([start])
            labels['extra'] = [n for n in self.label_names([c for c in chain if c != start])
                               if n not in names]
            strategy = {'output': {'detail': 'full'}, 'direction': arm, 'labels': labels,
                        'bounds': {'max_extension_bp': 110}}
            roll = rng.random()
            if roll < 0.5:
                strategy['branching'] = {'on_reconverge': 'keep'}
            elif roll < 0.8:
                strategy['exhaustive'] = True
            # else: merge, judged as merged arms are
            out.append({'id': self.next_id(), 'route': 'traverse', 'class': 'ample',
                        'request': {'seeds': [{'sequence': r.seq[p:p + L], 'labels': names}],
                                    'strategy': strategy}})
        return out


def generate(index, walker, config_name, seed, cases, n_random=40, n_caps=12):
    """The rule cases, the random requests, the cap probes, then the transition-chain requests
    (n_random // 5, at least 4; last, so the earlier ids are unchanged)."""
    g = TraverseGen(index, walker, seed, config_name, cases)
    return (g.rule_requests() + g.random_requests(n_random) + g.cap_requests(n_caps)
            + g.chain_requests(max(4, n_random // 5)))
