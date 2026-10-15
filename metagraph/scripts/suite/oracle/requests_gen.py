"""Deterministic requests (ORACLE-DESIGN.md section 5): /pattern requests over every axis of
phase 1, budget probes whose answers are legitimately incomplete, and /resolve requests.

Everything comes from one random.Random(seed): the same seed gives the same requests. Patterns
are drawn from the fixture records (windows at random and at the built-in features), their
reverse complements, IUPAC variants, translated frames and random strings.
"""

import random

from . import oracle as O

IUPAC_CODES = 'RYSWKMBDHVN'
IUPAC_COMP = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
              'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N'}
CODES = [1, 11, 1, 11, 2, 4, 6, 27, 31]


class Gen:
    def __init__(self, index, seed, config_name):
        self.index = index
        self.k = index.k
        self.rng = random.Random(f'{seed}:{config_name}')
        self.records = [r for r in index.records if len(r.seq) >= index.k]
        self.features = self._features()
        self.n = 0
        self.config = config_name

    # -- windows ---------------------------------------------------------------------------

    def _features(self):
        """(record, position) points the fixtures build in: SNP, element ends, repeats,
        homopolymer, planted substitutions, the N, palindromes, record ends."""
        pts = []
        by_name = {r.name: r for r in self.records}
        for r in self.records:
            pts += [(r, 0), (r, len(r.seq) - 1)]
            for motif in ('GAATTC', 'GGATCC', 'CCCGGG', 'ACGTACGT', 'ATATAT', 'GCGGCG',
                          'ACCGTTAACCGTTA', 'AAAAAAAAAAAA', 'N'):
                i = r.seq.find(motif)
                if i != -1:
                    pts.append((r, i + len(motif) // 2))
        if 'beta' in by_name:
            pts.append((by_name['beta'], 45))
        if 'eps' in by_name:
            pts += [(by_name['eps'], p) for p in (10, 40, 44, 70, 73, 77)]
        if 'alpha' in by_name:
            pts += [(by_name['alpha'], p) for p in (96, 146)]   # E's start and end in alpha
        return pts

    def window(self, L, feature=None):
        """A window of L bases of a record island (bases only), around a feature or anywhere."""
        return self.sourced_window(L, feature)[0]

    def sourced_window(self, L, feature=None):
        """(window, its record) as window() draws it; the record is None for a random string."""
        rng = self.rng
        for _ in range(200):
            if feature or rng.random() < 0.35:
                r, p = rng.choice(self.features)
                s = p - rng.randrange(0, L)
            else:
                r = rng.choice(self.records)
                s = rng.randrange(0, max(1, len(r.seq) - L + 1))
            s = max(0, min(s, len(r.seq) - L))
            w = r.seq[s:s + L]
            if len(w) == L and all(c in O.BASES for c in w):
                return w, r
        return ''.join(rng.choice(O.BASES) for _ in range(L)), None

    # -- patterns --------------------------------------------------------------------------

    def dna(self, L=None, rc=None):
        rng = self.rng
        L = L or (rng.choice([1, 2]) if rng.random() < 0.05 else rng.randint(3, self.k))
        roll = rng.random()
        if roll < 0.08:
            return ''.join(rng.choice(O.BASES) for _ in range(L))
        if roll < 0.16:
            h = ''.join(rng.choice(O.BASES) for _ in range(max(1, L // 2)))
            return h + O.revcomp(h)
        if roll < 0.22:
            r = rng.choice(self.records)
            seq = r.seq.split('N')[0] if rng.random() < 0.5 else r.seq.split('N')[-1]
            if len(seq) >= L:
                return seq[:L] if rng.random() < 0.5 else seq[-L:]
        w = self.window(L)
        if rc if rc is not None else rng.random() < 0.4:
            w = O.revcomp(w)
        if rng.random() < 0.05:
            w = w.lower()
        return w

    def iupac(self, L=None):
        rng = self.rng
        L = L or rng.randint(3, self.k)
        roll = rng.random()
        if roll < 0.15:
            # an even palindrome of codes (an odd centre would have to be its own complement)
            half = [rng.choice('ACGTRYSWKMBDHVN') for _ in range(max(1, L // 2))]
            p = half + [IUPAC_COMP[c] for c in reversed(half)]
            return ''.join(p) if p else 'WS'
        p = [c.upper() for c in self.dna(L)]
        L = len(p)          # a drawn palindrome is of even length
        for _ in range(rng.randint(1, min(3, L))):
            i = rng.randrange(L)
            if rng.random() < 0.6:
                p[i] = rng.choice([c for c in IUPAC_CODES if O.IUPAC[p[i]] <= O.IUPAC[c]])
            else:
                p[i] = rng.choice(IUPAC_CODES)
        if roll < 0.35:
            # an N run: leading (skipped by the search, and under a source dummy's '$' on a graph
            # without its mask), inner or trailing (the reverse orientation's leading run)
            run = rng.randint(1, max(1, L // 3))
            where = rng.choice(['lead', 'inner', 'trail'])
            if where == 'lead':
                p[:run] = ['N'] * run
            elif where == 'trail':
                p[L - run:] = ['N'] * run
            else:
                s = rng.randrange(0, max(1, L - run))
                p[s:s + run] = ['N'] * run
        return ''.join(p)

    def protein(self, code, max_res=None):
        rng = self.rng
        max_res = max_res or self.k // 3
        m = rng.randint(1, max_res)
        roll = rng.random()
        if roll < 0.1:
            return ''.join(rng.choice(O.AMINO_ACIDS + 'X') for _ in range(m))
        w = self.window(3 * m)
        if rng.random() < 0.4:
            w = O.revcomp(w)
        table = O.codon_table(code)
        res = [table[w[i:i + 3]] for i in range(0, 3 * m, 3)]
        for i in range(m):
            if res[i] == '*':
                continue
            x = rng.random()
            if x < 0.12:
                res[i] = 'X'
            elif x < 0.2:
                for amb, pair in O.AMBIGUOUS.items():
                    if res[i] in pair:
                        res[i] = amb
        return ''.join(res)

    def pattern(self, kind, code=1):
        if kind == 'dna':
            return {'dna': self.dna()}
        if kind == 'iupac':
            return {'iupac': self.iupac()}
        return {'protein': self.protein(code)}

    # -- requests --------------------------------------------------------------------------

    def next_id(self, prefix):
        self.n += 1
        return f'{prefix}-{self.config}-{self.n:04d}'

    def pattern_requests(self, count, labels_ok):
        rng = self.rng
        out = []
        for _ in range(count):
            mode = rng.choice(['count', 'all_or_count', 'all_or_count', 'partial'])
            req = {'mode': mode}
            scope = rng.choice(['any_offset', 'any_offset', 'suffix'])
            if scope != 'any_offset' or rng.random() < 0.3:
                req['scope'] = scope
            strands = rng.choice(['both', 'both', 'forward', 'reverse'])
            if strands != 'both' or rng.random() < 0.3:
                req['strands'] = strands
            m = rng.choice([0, 0, 0, 1, 1, 2, 3])
            code = rng.choice(CODES)
            kinds = [rng.choice(['dna', 'dna', 'iupac', 'protein']) for _ in range(rng.randint(1, 8))]
            pats = []
            for kind in kinds:
                p = self.pattern(kind, code)
                text = p[kind]
                L = 3 * len(text) if kind == 'protein' else len(text)
                if m and L < 2 * m + 2:
                    # a short pattern within m substitutions matches nearly every window
                    p = self.pattern('dna', code) if kind != 'protein' else p
                pats.append(p)
            if rng.random() < 0.15:
                pats.append(rng.choice([{'dna': ''}, {'dna': 'ACGU'}, {'iupac': 'ACXT'},
                                        {'protein': 'MOK'}]))
            if code in (27, 31) and rng.random() < 0.5:
                pats.append({'protein': 'M*'})      # '*' admits no codon there: no instance
            if any('protein' in p for p in pats) and (code != 1 or rng.random() < 0.5):
                req['genetic_code'] = code
            for i, p in enumerate(pats):
                if rng.random() < 0.3:
                    p['id'] = f'p{i}'
            req['patterns'] = [self.bounded(p, req, mode, m, code) for p in pats[:16]]
            if m:
                req['max_mismatches'] = m
            labels = labels_ok and not m and mode != 'count' and rng.random() < 0.5
            if labels_ok and not m and mode == 'count' and rng.random() < 0.15:
                labels = True
            if labels:
                req['output'] = {'labels': 'all'}
                if self.index.mode != 'basic':
                    req['allow_unbudgeted_annotation'] = True
            elif rng.random() < 0.3:
                req['output'] = {'labels': 'none'}
            if mode != 'count' and rng.random() < 0.2:
                req['max_contexts'] = rng.choice([1, 2, 3, 5, 8, 20])
            out.append({'id': self.next_id('P'), 'route': 'pattern', 'class': 'ample',
                        'request': req})
        return out

    MAX_RETRIEVED = 4000

    def bounded(self, p, req, mode, m, code):
        """|p|, or a redrawn pattern when a retrieval would list more than MAX_RETRIEVED
        contexts (a short degenerate pattern within m substitutions matches nearly every
        window: its list says little and its answer is megabytes)."""
        if mode == 'count':
            return p
        for _ in range(20):
            kind = next(k for k in ('dna', 'iupac', 'protein') if k in p)
            pat = O.Pattern(kind, p[kind], code)
            if pat.slot or pat.length > self.k:
                return p
            n = len(O.contexts(self.index, pat, req.get('scope', 'any_offset'),
                               req.get('strands', 'both'), m))
            if n <= self.MAX_RETRIEVED:
                return p
            q = self.pattern(kind if kind != 'iupac' else 'dna', code)
            q.update({k: v for k, v in p.items() if k == 'id'})
            p = q
        return p

    def budget_requests(self, count):
        """Requests a tight budget stops: legitimately incomplete answers (section 5.3)."""
        rng = self.rng
        out = []
        for i in range(count):
            mode = rng.choice(['count', 'all_or_count', 'partial'])
            req = {'mode': mode, 'patterns': []}
            if i % 3 == 2:
                # stop_at_threshold: a pattern with many contexts and a small threshold
                req['stop_at_threshold'] = True
                req['max_contexts'] = rng.choice([1, 2, 3, 5])
                req['patterns'] = [{'iupac': rng.choice(['NNA', 'ANNC', 'GNNNT', 'NRY'])},
                                   {'dna': self.dna(rng.randint(4, 7))}]
                if mode == 'count':
                    req['mode'] = 'partial'
            else:
                req['max_steps'] = rng.choice([1, 3, 10, 40, 150, 600])
                req['patterns'] = [{'iupac': self.iupac(rng.randint(5, self.k))}
                                   for _ in range(rng.randint(2, 4))]
                if rng.random() < 0.5:
                    req['max_mismatches'] = rng.choice([1, 2])
            out.append({'id': self.next_id('B'), 'route': 'pattern', 'class': 'budget',
                        'request': req})
        return out

    def periodic_runs(self):
        """(record, start, end, period) of the maximal periodic runs (period 1-7) long enough
        that whole periods can be deleted and k bases of the run remain."""
        k = self.k
        out = []
        for r in self.index.indexed:
            s = r.seq
            for p in range(1, 8):
                i = p
                while i < len(s):
                    j = i
                    while j < len(s) and s[j] == s[j - p] and s[j] in O.BASES:
                        j += 1
                    a = i - p
                    if j - a >= k + p and all(c in O.BASES for c in s[a:j]):
                        if not any(x[0] is r and x[1] <= a and j <= x[2] for x in out):
                            out.append((r, a, j, p))
                    i = j + 1
        return out

    def trace_requests(self, header_labels=True):
        """Queries that tell trace from k-mer presence (verification M2): a window across a
        periodic run with whole periods deleted from it. Every k-mer stays in the record, but
        the coordinate chain jumps where the query leaves the run, so a trace run breaks where
        a presence run does not."""
        rng = self.rng
        k = self.k
        out = []
        for r, a, e, p in self.periodic_runs():
            s = r.seq
            n = (e - a) // p
            j = rng.randint(1, max(1, (e - a - k) // p))          # periods deleted
            if (e - a) - j * p < k:
                continue
            fl = min(a, rng.randint(k, 20))
            fr = min(len(s) - e, rng.randint(k, 20))
            x = a + p * rng.randint(0, n - j)                     # where the deletion starts
            query = s[a - fl:x] + s[x + j * p:e + fr]
            if not all(c in O.BASES for c in query) or len(query) < k:
                continue
            labels = [r.column] + ([r.name] if header_labels else [])
            rng.shuffle(labels)
            out.append({'id': self.next_id('R'), 'route': 'resolve', 'class': 'ample',
                        'request': {'sequence': query, 'labels': labels, 'support': 'trace'}})
        return out

    def resolve_requests(self, count, header_labels=True):
        """Queries cut from the records (whole windows, junctions of two, a substituted base, a
        reverse complement, a window holding the record N or an inserted N); the labels mix the
        source records' column and header labels with others, in random order."""
        rng = self.rng
        k = self.k
        out = []
        columns = list(self.index.columns)
        headers = [r.name for r in self.index.indexed]
        for _ in range(count):
            roll = rng.random()
            sources = []
            if roll < 0.35:
                seq, r = self.sourced_window(rng.randint(k, 120))
                sources.append(r)
            elif roll < 0.55:
                a, ra = self.sourced_window(rng.randint(k, 60))
                b, rb = self.sourced_window(rng.randint(k, 60))
                seq = a + b
                sources += [ra, rb]
            elif roll < 0.7:
                w, r = self.sourced_window(rng.randint(2 * k, 120))
                s = list(w)
                i = rng.randrange(len(s))
                s[i] = rng.choice([b for b in O.BASES if b != s[i]])
                seq = ''.join(s)
                sources.append(r)
            elif roll < 0.8:
                w, r = self.sourced_window(rng.randint(k, 80))
                seq = O.revcomp(w)
                sources.append(r)
            elif roll < 0.9:
                r = next(r for r in self.index.records if 'N' in r.seq)
                i = r.seq.index('N')
                seq = r.seq[max(0, i - rng.randint(5, 40)):i + rng.randint(5, 40)]
                sources.append(r)
            else:
                w, r = self.sourced_window(rng.randint(2 * k, 80))
                i = rng.randrange(len(w))
                seq = w[:i] + 'N' + w[i + 1:]
                sources.append(r)
            labels = []
            for r in sources:
                if r is None:
                    continue
                labels.append(r.column)
                if header_labels and r.name in headers:
                    labels.append(r.name)
            pool = columns + (headers if header_labels else [])
            labels += rng.sample(pool, rng.randint(0, 3))
            if header_labels and rng.random() < 0.1:
                # a record shorter than k is not indexed: its header is refused (400)
                short = [r.name for r in self.index.records if len(r.seq) < k]
                if short:
                    labels.append(rng.choice(short))
            labels = list(dict.fromkeys(labels)) or [rng.choice(columns)]
            rng.shuffle(labels)
            req = {'sequence': seq, 'labels': labels}
            if header_labels and rng.random() < 0.3 and self.index.mode == 'basic':
                req['support'] = 'trace'
            if rng.random() < 0.2:
                req['min_block_kmers'] = rng.choice([2, 5, 20])
            out.append({'id': self.next_id('R'), 'route': 'resolve', 'class': 'ample',
                        'request': req})
        return out


def generate(index, config_name, seed, labels_ok, resolve_ok, header_labels, n_pattern=60,
             n_budget=12, n_resolve=30):
    g = Gen(index, seed, config_name)
    reqs = g.pattern_requests(n_pattern, labels_ok)
    reqs += g.budget_requests(n_budget)
    if resolve_ok:
        reqs += g.resolve_requests(n_resolve, header_labels)
        if n_resolve and index.mode == 'basic':
            reqs += g.trace_requests(header_labels)
    return reqs
