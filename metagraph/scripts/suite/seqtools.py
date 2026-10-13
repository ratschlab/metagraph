"""Pure sequence helpers of the suite (stdlib only): reverse complement, IUPAC sets, NCBI genetic codes,
translation and the information bits the host's floor gates on (SPEC-pattern-search.md §2, §12.2, §21.5).

Used by the builder (to make windows, variants and the expectations) and by checks.py (to recompute what an
answer must hold); nothing here touches the network or the disk.
"""

from __future__ import annotations

import math
import random
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

DNA = "ACGT"
_COMP = str.maketrans("ACGTRYSWKMBDHVNacgtryswkmbdhvn", "TGCAYRSWMKVHDBNtgcayrswmkvhdbn")

IUPAC: Dict[str, str] = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}
# the code of each base set, for making an IUPAC position that still holds a given base
CODE_OF_SET: Dict[str, str] = {"".join(sorted(v)): k for k, v in IUPAC.items()}


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def is_acgt(s: str) -> bool:
    return bool(s) and all(c in "ACGT" for c in s)


def iupac_matches(code: str, base: str) -> bool:
    return base in IUPAC[code]


def iupac_match(pattern: str, instance: str) -> bool:
    """Every base of |instance| lies in the set of the IUPAC code at its position."""
    return len(pattern) == len(instance) and all(b in IUPAC[c] for c, b in zip(pattern, instance))


def hamming_iupac(pattern: str, instance: str) -> List[int]:
    """The positions where |instance| lies outside |pattern|'s set (substitutions, §21.3)."""
    assert len(pattern) == len(instance)
    return [i for i, (c, b) in enumerate(zip(pattern, instance)) if b not in IUPAC[c]]


# ---------------------------------------------------------------------------------------------- genetic codes
# NCBI's genetic codes (gc.prt): ncbieaa = the residue of each of the 64 codons in NCBI's order (T, C, A, G; first
# base slowest), "*" = stop. Copied from the search service's pattern validation (pv_core.py at 20343460), which
# took them from taxdump gc.prt of 2025-09-20. Tables 27, 28 and 31 list the sense reading of their context stops,
# which is how the host reads them (SPEC-pattern-search.md §12.2, decision #21).
NCBIEAA: Dict[int, str] = {
    1: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    2: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
    3: "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    4: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    5: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
    6: "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    9: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    10: "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    11: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    12: "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    13: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
    14: "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    15: "FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    16: "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    21: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    22: "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    23: "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    24: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
    25: "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    26: "FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    27: "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    28: "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    29: "FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    30: "FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    31: "FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    32: "FFLLSSSSYY*WCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    33: "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
}
CODONS = [a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"]
AMINO = "ACDEFGHIKLMNPQRSTVWY"
# the ambiguity residues of the host (§12.2): X any sense codon, B = D|N, Z = E|Q, J = I|L
AMBIGUOUS = {"B": "DN", "Z": "EQ", "J": "IL"}


def codon_table(table: int) -> Dict[str, str]:
    return dict(zip(CODONS, NCBIEAA[int(table)]))


def translate(dna: str, table: int) -> str:
    m = codon_table(table)
    dna = dna.upper()
    return "".join(m.get(dna[i:i + 3], "X") for i in range(0, len(dna) - len(dna) % 3, 3))


def codons_of(residue: str, table: int) -> List[str]:
    """The codons the host admits for one residue in |table| (§12.2): X every sense codon, B/Z/J the union of
    their two residues, * the table's stop codons (none in 27, 28, 31)."""
    m = codon_table(table)
    residue = residue.upper()
    if residue == "X":
        return [c for c in CODONS if m[c] != "*"]
    if residue in AMBIGUOUS:
        return [c for c in CODONS if m[c] in AMBIGUOUS[residue]]
    return [c for c in CODONS if m[c] == residue]


def sense_substitutes(residue: str, table: int) -> List[str]:
    """The sense codons that code another residue (the host's substitutes of a residue, §21.3)."""
    own = set(codons_of(residue, table))
    m = codon_table(table)
    return [c for c in CODONS if m[c] != "*" and c not in own]


def has_stop(table: int) -> bool:
    return "*" in NCBIEAA[int(table)]


# ---------------------------------------------------------------------------------------------- information bits
def information_bits_iupac(pattern: str) -> float:
    """Σ log2(4 / |set_i|) (§2)."""
    return sum(math.log2(4 / len(IUPAC[c])) for c in pattern)


def information_bits_peptide(peptide: str, table: int) -> float:
    """6m − log2 of the codon strings the peptide admits (§12.2); a residue without codons counts as one exact
    codon (6 bits), as the host keeps the bits finite."""
    bits = 0.0
    for r in peptide:
        n = len(codons_of(r, table))
        bits += 6.0 if n == 0 else math.log2(64 / n)
    return bits


def mismatch_ball_size_iupac(pattern: str, m: int) -> int:
    """|B|: the strings over the pattern's positions within m substitutions (per position |set| matching and
    4 − |set| substituting bases; an N cannot substitute), §21.5."""
    # dp[t] = number of strings over the positions so far with exactly t substitutions
    dp = [1] + [0] * m
    for c in pattern:
        s = len(IUPAC[c])
        sub = 4 - s
        new = [0] * (m + 1)
        for t in range(m + 1):
            if dp[t] == 0:
                continue
            new[t] += dp[t] * s
            if t + 1 <= m and sub:
                new[t + 1] += dp[t] * sub
        dp = new
    return sum(dp)


def mismatch_information_bits_iupac(pattern: str, m: int) -> float:
    return 2 * len(pattern) - math.log2(mismatch_ball_size_iupac(pattern, m))


# ---------------------------------------------------------------------------------------------- generators
def random_dna(rng: random.Random, n: int, gc: float = 0.5) -> str:
    out = []
    for _ in range(n):
        if rng.random() < gc:
            out.append(rng.choice("GC"))
        else:
            out.append(rng.choice("AT"))
    return "".join(out)


def shuffled(rng: random.Random, s: str, avoid: Iterable[str] = ()) -> str:
    """A composition-preserving shuffle that differs from |s| and from every string in |avoid|."""
    avoid = set(avoid) | {s}
    chars = list(s)
    for _ in range(100):
        rng.shuffle(chars)
        t = "".join(chars)
        if t not in avoid:
            return t
    return "".join(chars)


def substitute_bases(window: str, positions: Sequence[int], rng: random.Random) -> Tuple[str, List[dict]]:
    """Another base at each position (seeded); returns the variant and the edits (pos, from, to)."""
    w = list(window)
    edits = []
    for p in sorted(positions):
        old = w[p]
        new = rng.choice([b for b in DNA if b != old])
        w[p] = new
        edits.append({"pos": p, "from": old, "to": new})
    return "".join(w), edits


def iupac_code_holding(base: str, size: int, rng: random.Random) -> str:
    """An IUPAC code of |size| bases (2, 3 or 4) whose set holds |base|."""
    choices = [k for k, v in IUPAC.items() if len(v) == size and base in v]
    return rng.choice(sorted(choices))


def spread_positions(rng: random.Random, n: int, t: int, avoid: Iterable[int] = ()) -> List[int]:
    """t distinct positions in [0, n), at least 2 apart where possible (spread substitutions)."""
    avoid = set(avoid)
    for _ in range(200):
        pos = sorted(rng.sample([i for i in range(n) if i not in avoid], t))
        if all(b - a >= 2 for a, b in zip(pos, pos[1:])):
            return pos
    return sorted(rng.sample([i for i in range(n) if i not in avoid], t))


def clustered_positions(rng: random.Random, n: int, t: int) -> List[int]:
    """t adjacent positions starting at a seeded offset (clustered substitutions)."""
    start = rng.randrange(0, n - t + 1)
    return list(range(start, start + t))


def low_complexity_dust_score(s: str) -> float:
    """The DUST score of a whole string: Σ c_t (c_t − 1) / 2 / (l − 1) over its l triplets (sdust's window score,
    without the perfect-interval search). Used by the builder only to keep a 'unique' window away from obvious
    low complexity; never as an expectation of the host's note."""
    if len(s) < 4:
        return 0.0
    counts: Dict[str, int] = {}
    for i in range(len(s) - 2):
        counts[s[i:i + 3]] = counts.get(s[i:i + 3], 0) + 1
    l = len(s) - 2
    return sum(c * (c - 1) / 2 for c in counts.values()) / (l - 1)


def find_all(haystack: str, needle: str) -> List[int]:
    out = []
    i = haystack.find(needle)
    while i != -1:
        out.append(i)
        i = haystack.find(needle, i + 1)
    return out


def occurrences_both_strands(record: str, s: str) -> List[Tuple[int, str]]:
    """(0-based start, strand) of |s| on the record: '+' where s is on the forward strand, '-' where rc(s) is."""
    out = [(p, "+") for p in find_all(record, s)]
    rc = revcomp(s)
    if rc != s:
        out += [(p, "-") for p in find_all(record, rc)]
    return sorted(out)


def first_codon_table_mismatch(dna: str, peptide: str, table: int) -> Optional[int]:
    t = translate(dna, table)
    for i, (a, b) in enumerate(zip(t, peptide)):
        if a != b:
            return i
    return None


# ---------------------------------------------------------------------------------------------- sdust
# A port of sdust (metagraph/external-libraries/sdust/sdust.c, minimap2's): the low-complexity intervals of a DNA
# string. The host states `low_complexity_pattern` on an exact pattern that sdust flags anywhere with the seeder's
# parameters T = 20, W = 64 (pattern_search.cpp is_low_complexity; SPEC-pattern-search §7.8), when its answer
# states no budget stop; the builder uses sdust_flags for that expectation (tests/test_seqtools.py checks the port
# against every stored /pattern answer of an exact DNA pattern).
_NT4 = {"A": 0, "C": 1, "G": 2, "T": 3, "a": 0, "c": 1, "g": 2, "t": 3}
_WLEN = 3
_WMSK = (1 << (_WLEN << 1)) - 1


def sdust(seq: str, T: int = 20, W: int = 64) -> List[Tuple[int, int]]:
    """The masked intervals [(start, finish)] of |seq|, as sdust_core computes them."""
    st = {"rv": 0, "rw": 0, "L": 0}
    cv, cw = [0] * 64, [0] * 64
    w = []
    P = []      # [start, finish, r, l], by descending start, then ascending finish
    res = []

    def save(start):
        if not P or P[-1][0] >= start:
            return
        p = P[-1]
        saved = False
        if res:
            s, f = res[-1]
            if p[0] <= f:
                saved = True
                res[-1] = (s, f if f > p[1] else p[1])
        if not saved:
            res.append((p[0], p[1]))
        i = len(P) - 1
        while i >= 0 and P[i][0] < start:
            i -= 1
        del P[i + 1:]

    def shift(t):
        if len(w) >= W - _WLEN + 1:
            s = w.pop(0)
            cw[s] -= 1
            st["rw"] -= cw[s]
            if st["L"] > len(w):
                st["L"] -= 1
                cv[s] -= 1
                st["rv"] -= cv[s]
        w.append(t)
        st["L"] += 1
        st["rw"] += cw[t]
        cw[t] += 1
        st["rv"] += cv[t]
        cv[t] += 1
        if cv[t] * 10 > T << 1:
            while True:
                s = w[len(w) - st["L"]]
                cv[s] -= 1
                st["rv"] -= cv[s]
                st["L"] -= 1
                if s == t:
                    break

    def find_perfect(start):
        c = cv[:]
        r = st["rv"]
        max_r = max_l = 0
        for i in range(len(w) - st["L"] - 1, -1, -1):
            t = w[i]
            r += c[t]
            c[t] += 1
            new_r, new_l = r, len(w) - i - 1
            if new_r * 10 > T * new_l:
                j = 0
                while j < len(P) and P[j][0] >= i + start:
                    p = P[j]
                    if max_r == 0 or p[2] * max_l > max_r * p[3]:
                        max_r, max_l = p[2], p[3]
                    j += 1
                if max_r == 0 or new_r * max_l >= max_r * new_l:
                    max_r, max_l = new_r, new_l
                    P.insert(j, [i + start, len(w) + (_WLEN - 1) + start, new_r, new_l])

    l = t = 0
    n = len(seq)
    for i in range(n + 1):
        b = _NT4.get(seq[i], 4) if i < n else 4
        if b < 4:
            l += 1
            t = ((t << 2) | b) & _WMSK
            if l >= _WLEN:
                start = (l - W if l - W > 0 else 0) + (i + 1 - l)
                save(start)
                shift(t)
                if st["rw"] * 10 > st["L"] * T:
                    find_perfect(start)
        else:
            start = (l - W + 1 if l - W + 1 > 0 else 0) + (i + 1 - l)
            while P:
                save(start)
                start += 1
            l = t = 0
    return res


def sdust_flags(seq: str, T: int = 20, W: int = 64) -> bool:
    """Whether sdust flags |seq| anywhere (the host reads it in overlapping pieces, which flag exactly when the
    whole does)."""
    return len(sdust(seq, T, W)) > 0
