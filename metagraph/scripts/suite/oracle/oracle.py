"""The brute-force oracle: the complete answer of a request, by plain enumeration over the records.

Deliberately simple and independent (ORACLE-DESIGN.md section 3): every k-mer of the index, every
offset of the scope and every searched orientation is tested position by position against the
pattern -- no interval pruning, no range search, no codon automaton, no regex, no pigeonhole
filter, no stopping logic. The IUPAC table, the reverse complement and the NCBI genetic codes are
this module's own, written from the SPEC (SPEC-pattern-search.md sections 2, 7.1-7.4, 12.2, 14,
21; SPEC-labeled-traversal-core.md section 4.2). It imports nothing of the engine, of checks.py,
seqtools.py, truth.py or the integration tests.
"""

import math

BASES = 'ACGT'
COMP = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

# the 15 IUPAC codes as base sets (SPEC section 4.2)
IUPAC = {
    'A': frozenset('A'), 'C': frozenset('C'), 'G': frozenset('G'), 'T': frozenset('T'),
    'R': frozenset('AG'), 'Y': frozenset('CT'), 'S': frozenset('CG'), 'W': frozenset('AT'),
    'K': frozenset('GT'), 'M': frozenset('AC'), 'B': frozenset('CGT'), 'D': frozenset('AGT'),
    'H': frozenset('ACT'), 'V': frozenset('ACG'), 'N': frozenset('ACGT'),
}

# NCBI translation tables (gc.prt 4.6), the ncbieaa column: the residue of each codon with the
# codons in TCAG order (TTT, TTC, TTA, TTG, TCT, ...). Written out here, never read from the
# engine.
GENETIC_CODES = {
    1: 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    2: 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG',
    4: 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    6: 'FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    11: 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    # no unconditional stop codon: TAA/TAG (and TGA) code a residue (SPEC section 12.2)
    27: 'FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    31: 'FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
}
CODONS = [a + b + c for a in 'TCAG' for b in 'TCAG' for c in 'TCAG']
AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'
PROTEIN_LETTERS = AMINO_ACIDS + 'XBZJ*'
AMBIGUOUS = {'B': 'DN', 'Z': 'EQ', 'J': 'IL'}

STRAND_NAMES = {'basic': {'fwd': '+', 'rev': '-', 'pal': '='},
                'canonical': {'fwd': 'forward', 'rev': 'reverse', 'pal': 'palindromic'},
                'primary': {'fwd': 'forward', 'rev': 'reverse', 'pal': 'palindromic'}}
# the order of orientations in results and in the entry's `strands` (SPEC section 7.9)
ORIENTATION_RANK = {'+': 0, 'forward': 0, '-': 1, 'reverse': 1, '=': 2, 'palindromic': 2}


def revcomp(s):
    return ''.join(COMP[c] for c in reversed(s))


def complement_set(bases):
    return frozenset(COMP[b] for b in bases)


def codon_table(code):
    return dict(zip(CODONS, GENETIC_CODES[code]))


def residue_codons(residue, code):
    """The codons a residue of a peptide admits in table |code| (SPEC section 12.2): an amino
    acid's own codons; X every sense codon; B, Z, J the codons of their two residues; '*' the
    table's stop codons (none in tables 27 and 31)."""
    table = codon_table(code)
    if residue == '*':
        return frozenset(c for c, aa in table.items() if aa == '*')
    if residue == 'X':
        return frozenset(c for c, aa in table.items() if aa != '*')
    allowed = AMBIGUOUS.get(residue, residue)
    return frozenset(c for c, aa in table.items() if aa in allowed and aa != '*')


def sense_codons(code):
    return frozenset(c for c, aa in codon_table(code).items() if aa != '*')


def colex_key(kmer):
    """The BOSS edge order of real k-mers: the source node (k - 1 bases) co-lexicographically,
    then the edge label (SPEC section 7.9, ORACLE-DESIGN.md section 3.3)."""
    return (kmer[:-1][::-1], kmer[-1])


class Index:
    """The k-mers of an index built from |records| at |k| in graph mode |mode|, with the label
    and placement facts of its annotation (one column per fixture column)."""

    def __init__(self, records, k, mode):
        assert mode in ('basic', 'canonical', 'primary')
        self.records = list(records)
        self.k = k
        self.mode = mode
        self.forward = set()
        for r in self.records:
            for island in _islands(r.seq):
                for s in range(len(island) - k + 1):
                    self.forward.add(island[s:s + k])
        if mode == 'basic':
            self.kmers = set(self.forward)
        else:
            self.kmers = self.forward | {revcomp(x) for x in self.forward}
        self.ordered = sorted(self.kmers, key=colex_key)
        self.rank = {x: i for i, x in enumerate(self.ordered)}
        self.names = STRAND_NAMES[mode]
        self.columns = []
        for r in self.records:
            if r.column not in self.columns:
                self.columns.append(r.column)
        # The record mapping (.seqs) lists a column's records of at least k bases, in FASTA
        # order: a record shorter than k is not indexed (capabilities records_shorter_than_k
        # "not_indexed"), has no header label and takes no seq_id. A record's k-mer coordinates
        # in its column follow the indexed records before it, len - k + 1 slots each (every
        # start, an N-holding one included). ORACLE-DESIGN.md section 3.4.
        self.seq_ids = {}
        self.column_offset = {}
        nxt, off = {}, {}
        for r in self.records:
            if len(r.seq) < k:
                continue
            self.seq_ids[(r.column, r.name)] = nxt.get(r.column, 0)
            nxt[r.column] = nxt.get(r.column, 0) + 1
            self.column_offset[(r.column, r.name)] = off.get(r.column, 0)
            off[r.column] = off.get(r.column, 0) + len(r.seq) - k + 1
        self.indexed = [r for r in self.records if len(r.seq) >= k]

    # -- labels and placement (SPEC section 14) -------------------------------------------

    def starts_in_record(self, kmer, record):
        """The 0-based starts of |kmer| on the forward strand of |record|."""
        out = []
        i = record.seq.find(kmer)
        while i != -1:
            out.append(i)
            i = record.seq.find(kmer, i + 1)
        return out

    def labels_of(self, kmer):
        """The columns annotating |kmer|'s row: on BASIC the records holding it on their forward
        strand; on CANONICAL and PRIMARY those holding it or its reverse complement (one row
        serves both orientations)."""
        cache = self.__dict__.setdefault('_labels', {})
        if kmer not in cache:
            cache[kmer] = self._labels_of(kmer)
        return cache[kmer]

    def _labels_of(self, kmer):
        out = set()
        for r in self.records:
            if kmer in r.seq or (self.mode != 'basic' and revcomp(kmer) in r.seq):
                out.add(r.column)
        return out

    def occurrences(self, kmer, column, offset, length, strand):
        """The placed occurrences of column |column| in the context (strand, kmer, offset) on a
        BASIC index with record placement: (seq_id, record, strand, nt_coords, nt_length)."""
        out = []
        for r in self.records:
            if r.column != column:
                continue
            for s in self.starts_in_record(kmer, r):
                start = s + offset + 1
                out.append((self.seq_ids[(r.column, r.name)], r.name, strand,
                            f'{start}-{start + length - 1}', len(r.seq)))
        return sorted(out, key=lambda o: (o[0], int(o[3].split('-')[0]), o[2]))

    # -- /resolve (SPEC-labeled-traversal-core.md section 4.2) ------------------------------

    def in_graph(self, kmer):
        if self.mode == 'basic':
            return kmer in self.forward
        return kmer in self.kmers

    def record_holds(self, record, kmer):
        if kmer in record.seq:
            return True
        return self.mode != 'basic' and revcomp(kmer) in record.seq


def _islands(seq):
    """The maximal runs of A, C, G, T of a record (a DNA4 build splits a record at every other
    symbol and keeps no k-mer across one)."""
    out, cur = [], []
    for c in seq.upper():
        if c in BASES:
            cur.append(c)
        else:
            if cur:
                out.append(''.join(cur))
            cur = []
    if cur:
        out.append(''.join(cur))
    return out


# -- patterns -----------------------------------------------------------------------------

class Pattern:
    """A request pattern as the oracle reads it: its kind, its text in upper case, its length L
    in bases, and the matcher of one oriented reading."""

    def __init__(self, kind, text, genetic_code=1):
        self.kind = kind
        self.raw = text
        self.text = text.upper()
        self.code = genetic_code
        self.slot = None
        if kind == 'dna':
            ok = self.text and all(c in BASES for c in self.text)
        elif kind == 'iupac':
            ok = self.text and all(c in IUPAC for c in self.text)
        else:
            ok = self.text and all(c in PROTEIN_LETTERS for c in self.text)
        if not ok:
            self.slot = 'bad_alphabet'
            return
        if kind == 'protein':
            self.residues = list(self.text)
            self.length = 3 * len(self.residues)
            self.sets = [residue_codons(r, genetic_code) for r in self.residues]
            self.sense = sense_codons(genetic_code)
            self.table = codon_table(genetic_code)
            self.no_instance = any(not s for s in self.sets)
        else:
            self.length = len(self.text)
            self.sets = [IUPAC[c] for c in self.text]
            self.no_instance = False

    @property
    def palindromic(self):
        if self.kind == 'protein':
            m = len(self.sets)
            return all(self.sets[i] == frozenset(revcomp(c) for c in self.sets[m - 1 - i])
                       for i in range(m))
        L = self.length
        return all(self.sets[i] == complement_set(self.sets[L - 1 - i]) for i in range(L))

    def substitutions(self, r, max_mismatches):
        """The substitutions of the reading |r| (len L, bases only) against the pattern, in the
        pattern's coordinates, or None when |r| is no instance within |max_mismatches|.
        DNA/IUPAC: (position, expected code, found base); a position whose set holds the base
        holds none (N never). Peptides (per residue, sense substitutes, SPEC section 21.3):
        (residue, expected residue, found residue, codon); X and '*' are never substituted, a
        stop codon never substitutes."""
        subs = []
        if self.kind == 'protein':
            for i, allowed in enumerate(self.sets):
                codon = r[3 * i:3 * i + 3]
                if codon in allowed:
                    continue
                res = self.residues[i]
                if res in ('X', '*') or codon not in self.sense:
                    return None
                subs.append((i, res, self.table[codon], codon))
                if len(subs) > max_mismatches:
                    return None
            return subs
        for i, allowed in enumerate(self.sets):
            if r[i] in allowed:
                continue
            subs.append((i, self.text[i], r[i]))
            if len(subs) > max_mismatches:
                return None
        return subs

    # -- information (SPEC sections 2, 12.2, 21.5) ------------------------------------------

    def _factors(self):
        """Per position (residue): (strings matching, strings substituting)."""
        if self.kind == 'protein':
            out = []
            for res, allowed in zip(self.residues, self.sets):
                if not allowed:
                    out.append((1, 0))            # counted as one exact codon (SPEC 12.2)
                elif res in ('X', '*'):
                    out.append((len(allowed), 0))
                else:
                    out.append((len(allowed), len(self.sense - allowed)))
            return out
        return [(len(s), 4 - len(s)) for s in self.sets]

    def information_bits(self):
        unit = 64 if self.kind == 'protein' else 4
        return sum(math.log2(unit / a) for a, _ in self._factors())

    def mismatch_information_bits(self, m):
        """2L - log2 |B|, B the strings within m substitutions (by levels, position by
        position)."""
        levels = [1] + [0] * m
        for a, b in self._factors():
            levels = [levels[t] * a + (levels[t - 1] * b if t else 0) for t in range(m + 1)]
        return 2 * self.length - math.log2(sum(levels))


# -- the contexts of a pattern ------------------------------------------------------------

class Context:
    __slots__ = ('orientation', 'kmer', 'offset', 'level', 'subs')

    def __init__(self, orientation, kmer, offset, level, subs):
        self.orientation = orientation
        self.kmer = kmer
        self.offset = offset
        self.level = level
        self.subs = subs

    @property
    def identity(self):
        return (self.orientation, self.kmer, self.offset)

    def __repr__(self):
        return f'Context({self.orientation}, {self.kmer}, {self.offset}, level {self.level})'


def orientations(index, pattern, strands):
    """The oriented readings searched (SPEC section 7.3): a palindrome once, whatever |strands|
    says; else forward and/or reverse."""
    names = index.names
    if pattern.palindromic:
        return [names['pal']]
    out = []
    if strands in ('both', 'forward'):
        out.append(names['fwd'])
    if strands in ('both', 'reverse'):
        out.append(names['rev'])
    return out


def offsets(k, L, scope):
    return [k - L] if scope == 'suffix' else list(range(0, k - L + 1))


def contexts(index, pattern, scope='any_offset', strands='both', max_mismatches=0):
    """Every graph context (orientation, k-mer, offset) of |pattern| (L <= k) in the scope and
    strands, each with its level and substitutions, in the answer order of SPEC sections 7.9 and
    21.7 on a DBGSuccinct (level, BOSS order, offset, orientation)."""
    k, L = index.k, pattern.length
    assert L <= k, 'patterns longer than k are out of phase 1'
    out = []
    if pattern.no_instance:
        return out
    oris = orientations(index, pattern, strands)
    offs = offsets(k, L, scope)
    for kmer in index.ordered:
        for p in offs:
            w = kmer[p:p + L]
            for o in oris:
                r = revcomp(w) if o in ('-', 'reverse') else w
                subs = pattern.substitutions(r, max_mismatches)
                if subs is not None:
                    out.append(Context(o, kmer, p, len(subs), subs))
    out.sort(key=lambda c: (c.level, index.rank[c.kmer], c.offset, ORIENTATION_RANK[c.orientation]))
    return out


def counts(index, pattern, ctxs, scope, strands, max_mismatches):
    """The exact counts of a complete answer: total, suffix, by_offset, by_strand or
    by_orientation, by_mismatches (None when max_mismatches is 0)."""
    k, L = index.k, pattern.length
    by_offset = {str(p): 0 for p in offsets(k, L, scope)}
    by_orient = {}
    for o in orientations(index, pattern, strands):
        key = 'both' if o == '=' else o
        by_orient[key] = 0
    by_mm = {str(t): 0 for t in range(max_mismatches + 1)} if max_mismatches else None
    for c in ctxs:
        by_offset[str(c.offset)] += 1
        by_orient['both' if c.orientation == '=' else c.orientation] += 1
        if by_mm is not None:
            by_mm[str(c.level)] += 1
    return {'total': len(ctxs), 'suffix': sum(1 for c in ctxs if c.offset == k - L),
            'by_offset': by_offset, 'by_orient': by_orient, 'by_mismatches': by_mm}


# -- the upper bound U of a graph without its dummy-edge mask (SPEC section 7.4) ----------------

def source_dummies(index):
    """The source dummy k-mers of the BOSS graph of the index: '$' * j + x[:k - j], 1 <= j < k,
    for every node x (the first k - 1 bases of a k-mer) that no k-mer enters (no k-mer ends in
    x). The sink dummies (W = '$') and '$' * k never hold a pattern base."""
    cache = index.__dict__.get('_dummies')
    if cache is None:
        k = index.k
        nodes = {x[:k - 1] for x in index.kmers}
        entered = {x[1:] for x in index.kmers}
        cache = sorted({'$' * j + x[:k - j] for x in nodes - entered for j in range(1, k)})
        index._dummies = cache
    return cache


def leading_n_run(pattern, orientation):
    """The leading N run of the oriented pattern that the search skips on a $ACGT graph (SPEC
    section 7.8): P's for forward and palindromic, P's trailing run for reverse; a window of N
    only keeps its last position; none for a peptide (X excludes the stops)."""
    if pattern.kind == 'protein':
        return 0
    text = pattern.text if orientation not in ('-', 'reverse') else pattern.text[::-1]
    n = len(text) - len(text.lstrip('N'))
    return min(n, pattern.length - 1)


def dummy_contexts(index, pattern, scope='any_offset', strands='both', max_mismatches=0):
    """The contexts a graph without its mask counts in source dummies besides its k-mers
    (U - the masked count): a dummy holds the oriented pattern at an offset where its '$' run
    ends at or before the pattern's skipped leading N run, the rest matching within
    |max_mismatches| ('$' is never a substitute). Contexts of the dummy strings, each once."""
    k, L = index.k, pattern.length
    out = []
    if pattern.no_instance:
        return out
    oris = orientations(index, pattern, strands)
    for d in source_dummies(index):
        j = len(d) - len(d.lstrip('$'))
        for p in offsets(k, L, scope):
            w = d[p:p + L]
            dollars = max(0, j - p)
            for o in oris:
                if dollars > leading_n_run(pattern, o):
                    continue
                if dollars == 0:
                    r = revcomp(w) if o in ('-', 'reverse') else w
                    subs = pattern.substitutions(r, max_mismatches)
                else:
                    # the '$' lie under the leading N run: the rest is read in the orientation
                    # (for reverse, against the complemented sets of P read backwards)
                    subs = _oriented_substitutions(pattern, o, w, dollars, max_mismatches)
                if subs is not None:
                    out.append(Context(o, d, p, len(subs), subs))
    return out


def _oriented_substitutions(pattern, orientation, w, dollars, m):
    L = pattern.length
    subs = []
    for i in range(dollars, L):
        if orientation in ('-', 'reverse'):
            allowed = complement_set(pattern.sets[L - 1 - i])
        else:
            allowed = pattern.sets[i]
        if w[i] not in allowed:
            subs.append(i)
            if len(subs) > m:
                return None
    return subs


# -- labels over a list of returned contexts ------------------------------------------------

def label_view(index, pattern, returned):
    """The label projection over the returned contexts (SPEC sections 14.2-14.6): per context
    its labels in label order with their occurrences, by_label, counts.labels and
    counts.occurrences. Occurrences only on BASIC (record placement)."""
    L = pattern.length
    k = index.k
    record = index.mode == 'basic'
    per_ctx = []
    column_ctx = {}
    column_suffix = {}
    column_occ = {}
    for c in returned:
        cols = index.labels_of(c.kmer)
        occ = {}
        for col in cols:
            column_ctx[col] = column_ctx.get(col, 0) + 1
            if c.offset == k - L:
                column_suffix[col] = column_suffix.get(col, 0) + 1
            if record:
                o = index.occurrences(c.kmer, col, c.offset, L, c.orientation)
                occ[col] = o
                column_occ.setdefault(col, set()).update((x[0], x[3], x[2]) for x in o)
        per_ctx.append((c, cols, occ))
    order = sorted(column_ctx, key=lambda col: (-column_ctx[col], col))
    rank = {col: i for i, col in enumerate(order)}
    results = []
    for c, cols, occ in per_ctx:
        labels = sorted(cols, key=lambda col: rank[col])
        results.append({'labels': labels, 'occurrences': {col: occ.get(col) for col in labels}})
    by_label = [{'column': col, 'contexts': column_ctx[col],
                 'contexts_suffix': column_suffix.get(col, 0),
                 'occurrences': len(column_occ.get(col, ())) if record else None}
                for col in order]
    total_occ = sum(len(v) for v in column_occ.values()) if record else None
    return {'results': results, 'by_label': by_label, 'labels': len(order),
            'occurrences': total_occ}


# -- /resolve -----------------------------------------------------------------------------

def runs_of(flags):
    """The maximal runs [a, b) of True in |flags|."""
    out, start = [], None
    for i, f in enumerate(flags):
        if f and start is None:
            start = i
        elif not f and start is not None:
            out.append([start, i])
            start = None
    if start is not None:
        out.append([start, len(flags)])
    return out


def resolve(index, sequence, labels, min_block_kmers=1, support='kmer'):
    """The profile of explicit labels over |sequence| (SPEC-labeled-traversal-core.md 4.2):
    num_kmers, graph_runs, per label (in request order) kind, kmers_supported and runs, and the
    candidates. A label is a column name or a record header (kind 'column' / 'header'); an
    unknown name -> ValueError."""
    k = index.k
    seq = sequence.upper()
    n = max(0, len(seq) - k + 1)
    kmers = [seq[i:i + k] for i in range(n)]
    present = [all(c in BASES for c in x) and index.in_graph(x) for x in kmers]
    columns = set(index.columns)
    headers = {r.name: r for r in index.indexed}
    out_labels = []
    for name in labels:
        if name in columns:
            kind = 'column'
            recs = [r for r in index.records if r.column == name]
        elif name in headers:
            kind = 'header'
            recs = [headers[name]]
        else:
            raise ValueError(f'unknown label {name!r}')
        if support == 'trace':
            runs = _trace_runs(index, kmers, present, recs, kind)
            supported = sum(b - a for a, b in runs)
        else:
            flags = [present[i] and any(index.record_holds(r, kmers[i]) for r in recs)
                     for i in range(n)]
            runs = runs_of(flags)
            supported = sum(flags)
        out_labels.append({'label': name, 'kind': kind, 'kmers_supported': supported,
                           'runs': runs})
    groups = {}
    for pos, lab in enumerate(out_labels):
        for a, b in lab['runs']:
            if b - a >= min_block_kmers:
                groups.setdefault((a, b), []).append(lab['label'])
    candidates = [{'kmer_interval': [a, b], 'bp_interval': [a, b + k - 1], 'labels': names}
                  for (a, b), names in groups.items()]
    candidates.sort(key=lambda c: (-(c['kmer_interval'][1] - c['kmer_interval'][0]),
                                   -len(c['labels']), c['kmer_interval'][0]))
    return {'k': k, 'num_kmers': n, 'graph_runs': runs_of(present), 'labels': out_labels,
            'candidates': candidates}


def _trace_runs(index, kmers, present, records, kind):
    """Trace support (SPEC-labeled-traversal-core.md 4.2): the label's coordinates at k-mer i
    that continue the run's chain (c - 1 among the chain's at i - 1) carry it on; when none
    does, the run ends and a new one starts at i with all of i's coordinates. A header label's
    coordinates are its record's local starts; a column label's are column coordinates (the
    column's indexed records concatenated, so a chain can cross a record boundary: the stated
    limitation trace_record_boundaries)."""
    def coords(x):
        out = set()
        for r in records:
            if len(r.seq) < index.k:
                continue
            base = 0 if kind == 'header' else index.column_offset[(r.column, r.name)]
            out.update(base + s for s in index.starts_in_record(x, r))
        return out
    runs, start, live = [], None, set()
    for i, x in enumerate(kmers):
        pos = coords(x) if present[i] else set()
        if not pos:
            if start is not None:
                runs.append([start, i])
            start, live = None, set()
            continue
        cont = {p for p in pos if p - 1 in live}
        if start is not None and cont:
            live = cont
            continue
        if start is not None:
            runs.append([start, i])
        start, live = i, pos
    if start is not None:
        runs.append([start, len(kmers)])
    return runs
