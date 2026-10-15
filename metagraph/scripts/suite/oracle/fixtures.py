"""The fixture records: a few hand-designed tiny records, the same bytes on every run.

Each feature the oracle suite needs is built in by construction (ORACLE-DESIGN.md section 1):
a SNP bubble (branching and reconvergence), a shared element followed by divergent blocks
(branching), the reverse complement of a shared block, a period-7 tandem repeat (a cycle at
k = 11), a homopolymer (a self-loop k-mer), even palindromes, planted mismatch neighbours (one,
two and three substitutions inside one k-window), a record split by an N, records shorter than k
(last and first in their column), and columns of several records.
"""

import random

SEED = 20261015

# the order the columns are annotated in (their column ids): deliberately not their name order,
# so that a label order taken from column ids cannot pass as one taken from names
COLUMN_ORDER = ['c_zeta', 'c_alpha', 'c_theta', 'c_delta', 'c_beta', 'c_gamma', 'c_eps']

_COMP = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}


def _rc(s):
    return ''.join(_COMP[c] for c in reversed(s))


class Record:
    """One FASTA record of a column: its header, its text, its index in the column's file."""

    def __init__(self, column, seq_id, name, seq):
        self.column = column
        self.seq_id = seq_id
        self.name = name
        self.seq = seq

    def __repr__(self):
        return f'Record({self.column!r}, {self.seq_id}, {self.name!r}, {len(self.seq)} bp)'


def _substitute(seq, positions, rng):
    s = list(seq)
    for p in positions:
        s[p] = rng.choice([b for b in 'ACGT' if b != s[p]])
    return ''.join(s)


def build_records(seed=SEED):
    """{column: [Record]} in COLUMN_ORDER, each column's records in FASTA order."""
    rng = random.Random(seed)

    def block(n):
        return ''.join(rng.choice('ACGT') for _ in range(n))

    A, E, B, C = block(90), block(50), block(90), block(90)
    D, F, G, H = block(70), block(80), block(80), block(60)
    R = [block(n) for n in (30, 30, 20, 20, 10, 40, 100, 60)]
    unit = 'ACCGTTA'                      # period 7 < k: a cycle of the de Bruijn graph at k = 11
    # A with its middle base changed: beta and alpha share A's flanks, a bubble around the SNP
    a_snp = A[:45] + {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}[A[45]] + A[46:]

    alpha = A[:30] + 'GAATTC' + A[30:] + E + B
    beta = a_snp + E + C
    gamma = _rc(E) + unit * 6 + D + 'A' * 24 + H
    # planted neighbours of alpha[100:200]: an isolated substitution, a pair and a triple within
    # one k-window
    eps = _substitute(alpha[100:200], [10, 40, 44, 70, 73, 77], rng)
    # palindromes, among them a palindromic 10-mer and 12-mer (k-mers equal to their reverse
    # complement at an even k)
    zeta = ('ACGT' * 3 + R[0] + 'GGATCC' + R[1] + 'AT' * 8 + R[2] + 'GCG' * 4 + R[3]
            + 'CCCGGG' + R[4] + 'ACGTTAACGT' + R[2][:7] + 'GAATTCGAATTC')
    delta1 = F + 'N' + G                  # an N splits the record: no k-mer spans it
    delta2 = _rc(B[10:80]) + R[5]         # the reverse complement of a block alpha holds
    delta3 = 'ACGTTGCA'                   # shorter than k: no k-mer; the last of its column
    theta1 = 'GATTACA'                    # shorter than k, the FIRST of its column
    theta2 = R[6] + E[:20] + R[7]

    records = {
        'c_alpha': [('alpha', alpha)],
        'c_beta': [('beta', beta)],
        'c_gamma': [('gamma', gamma)],
        'c_eps': [('eps', eps)],
        'c_zeta': [('zeta', zeta)],
        'c_delta': [('delta1', delta1), ('delta2', delta2), ('delta3', delta3)],
        'c_theta': [('theta1', theta1), ('theta2', theta2)],
    }
    out = {}
    for column in COLUMN_ORDER:
        out[column] = [Record(column, i, name, seq) for i, (name, seq) in enumerate(records[column])]
    return out


def all_records(columns):
    """Every record, in column order then FASTA order."""
    return [r for column in COLUMN_ORDER if column in columns for r in columns[column]]


def write_fasta(columns, directory):
    """One FASTA file per column, <directory>/<column>.fa; -> {column: path}."""
    import os
    paths = {}
    for column, records in columns.items():
        path = os.path.join(directory, column + '.fa')
        with open(path, 'w') as f:
            for r in records:
                f.write(f'>{r.name}\n{r.seq}\n')
        paths[column] = path
    return paths
