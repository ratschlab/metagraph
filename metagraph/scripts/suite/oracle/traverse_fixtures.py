"""The traversal fixture (ORACLE-DESIGN.md section 10.1): records built so that every walker rule
meets one case it permits and one it forbids, with analytically known walks. The same bytes on every
run; one column per record (t_<record>) except t_amb, which holds two.

Cases (k = 11), each with its seed and labels in CASES:
- divergence: div1 = S1 + X1, div2 = S1 + X2 (two labels part at the seed boundary);
- ambiguous: column t_amb = {amb1 = S2 + Y1, amb2 = S2 + Y2} (one label on both sides);
- cycle: cyc = P3 + U x 5 + Q3, U of period 7;
- bubble: bub1 = S4 + M + T4, bub2 = S4 + M' + T4 (one SNP, reconverging k bases later);
- transition: sw1 = S5 + W1, sw2 = W1[20:] + W2, sw3 = W2[15:] + W4 (a chain of two switches);
- fork: fa = S7 + V1, fb = V1[15:] + V2, fc = V1[15:] + V3 (a switch into either of two successors);
- quorum: qa, qb, qc = S6 + Z1, qd = S6 + Z2 (three labels on one side, one on the other);
- circle: circ = P8 + Q8 + P8 (the walk re-enters the seed);
- homopolymer: hom = P9 + A x 14 + Q9 (a self-loop k-mer);
- hairpin: hp = P10 + GAATTCGAATTC + Q10 (a self-reverse-complementary 12-mer: a hairpin step at
  k = 11 on CANONICAL and PRIMARY);
- trace: tr = J1 + K1 + J1[-15:] + K2 (k-mer presence is ambiguous after J1, the coordinates are not);
- trace_break: tb = J2[-15:] + K5 + J2 (the seed J2 ends the record: its chain ends there while J2's
  tail goes on into K5 earlier in the record).
"""

import random

from .fixtures import Record

SEED = 20261016
K = 11


def _rc(s):
    comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join(comp[c] for c in reversed(s))


def build():
    """-> ({column: [Record]} in column order, {case: (seed, labels, extra)})."""
    rng = random.Random(SEED)

    def block(n, avoid_first=None):
        while True:
            s = ''.join(rng.choice('ACGT') for _ in range(n))
            if avoid_first is None or s[0] not in avoid_first:
                return s

    S1, S2, S4, S5, S6, S7 = (block(30) for _ in range(6))
    X1 = block(30)
    X2 = block(30, avoid_first=X1[0])
    Y1 = block(30)
    Y2 = block(30, avoid_first=Y1[0])
    P3, Q3 = block(30), block(30)
    U = 'TGCAGGA'
    M = block(21)
    M2 = M[:10] + {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}[M[10]] + M[11:]
    T4 = block(30)
    W1, W2, W4 = block(40), block(30), block(30)
    V1 = block(40)
    V2 = block(30)
    V3 = block(30, avoid_first=V2[0])
    Z1 = block(30)
    Z2 = block(30, avoid_first=Z1[0])
    P8, Q8 = block(20), block(30)
    P9, Q9 = block(30), block(30)
    P10, Q10 = block(30), block(30)
    J1 = block(30)
    K1 = block(20)
    K2 = block(20, avoid_first=K1[0])
    J2, K5 = block(30), block(20)
    records = [
        ('t_div1', 'div1', S1 + X1), ('t_div2', 'div2', S1 + X2),
        ('t_amb', 'amb1', S2 + Y1), ('t_amb', 'amb2', S2 + Y2),
        ('t_cyc', 'cyc', P3 + U * 5 + Q3),
        ('t_bub1', 'bub1', S4 + M + T4), ('t_bub2', 'bub2', S4 + M2 + T4),
        ('t_sw1', 'sw1', S5 + W1), ('t_sw2', 'sw2', W1[20:] + W2),
        ('t_sw3', 'sw3', W2[15:] + W4),
        ('t_fa', 'fa', S7 + V1), ('t_fb', 'fb', V1[15:] + V2), ('t_fc', 'fc', V1[15:] + V3),
        ('t_qa', 'qa', S6 + Z1), ('t_qb', 'qb', S6 + Z1), ('t_qc', 'qc', S6 + Z1),
        ('t_qd', 'qd', S6 + Z2),
        ('t_circ', 'circ', P8 + Q8 + P8),
        ('t_hom', 'hom', P9 + 'A' * 14 + Q9),
        ('t_hp', 'hp', P10 + 'GAATTCGAATTC' + Q10),
        ('t_tr', 'tr', J1 + K1 + J1[-15:] + K2),
        ('t_tb', 'tb', J2[-15:] + K5 + J2),
    ]
    # annotated in reverse name order: column ids are not the names' order
    order = sorted({c for c, _, _ in records}, reverse=True)
    columns = {c: [] for c in order}
    for c, name, seq in records:
        columns[c].append(Record(c, len(columns[c]), name, seq))
    cases = {
        'divergence': (S1, ['div1', 'div2'], []),
        'ambiguous': (S2, ['t_amb'], []),
        'cycle': (P3[-15:], ['cyc'], []),
        'bubble': (S4, ['bub1', 'bub2'], []),
        'transition': (S5, ['sw1'], ['sw2', 'sw3']),
        'fork': (S7, ['fa'], ['fb', 'fc']),
        'quorum': (S6, ['qa', 'qb', 'qc', 'qd'], []),
        'circle': (P8, ['circ'], []),
        'homopolymer': (P9[-15:], ['hom'], []),
        'hairpin': (P10[-15:], ['hp'], []),
        'trace': (J1, ['tr'], []),
        'trace_break': (J2, ['tb'], []),
    }
    return columns, cases


def column_order(columns):
    return list(columns)
