"""What the area generators share: the context (records, seeded randomness), window picking, the request
builders of the three routes (every default named explicitly, so that the direct and the MCP arm ask the same
question), the MCP mapping of DESIGN §4 with the service-by-design case ids of the async session, and the
scenario constructor of DESIGN §8.1.

Pure apart from reading the truth records it is given; no network.
"""

from __future__ import annotations

import copy
import hashlib
import json
import os
import random
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import seqtools as st
from truth import K, TRUTH_DIR, Record

DB = "refseq33m-experimental"

# Every pattern request names these (async finding 1: the MCP defaults differ: mode count, max_contexts 1,000)
PATTERN_TIME_BUDGET_MS = 60000
PATTERN_MAX_CONTEXTS = 10000
PATTERN_MAX_PREDICATE_CONTEXTS = 10000   # the MCP's default and the staging host's cap
# The host's defaults (the capabilities: default_scope any_offset, default_strands both, default_projection none,
# default_occurrences true; SPEC-pattern-search.md §4.1: predicate_strands "either"), equal to the MCP's
# (async's table), named in every body where they act (pattern_request)
PATTERN_DEFAULTS = {"scope": "any_offset", "strands": "both", "labels": "none", "occurrences": True,
                    "predicate_strands": "either"}
# Explicit-label resolves name the deadline and the query bound the service sends (its own max_query_bp 2,000 / 500
# for a discovery), so that both arms send the same body. The deadline is the host's cap, 60,000 ms on staging (the
# snapshot's resolve.time_budget.max_time_ms; the service's own default, TRAVERSAL_RESOLVE_TIME_BUDGET_MS, held to
# it): v4's 20,000 ms met cold first calls of 21-37 s under the Logan load (a relayed 503 or the service's timeout,
# triage-v4 G22, G23), which 60,000 ms holds with margin; a scenario that tests a budget names its own
RESOLVE_TIME_BUDGET_MS = 60000
RESOLVE_MAX_QUERY_BP = 2000
DISCOVER_MAX_QUERY_BP = 500
# The MCP's walk and trie presets (async's table, traverse.unnamed_strategy_fields; the budgets as the live
# traverse_capabilities states them on 2026-10-13: walk 15,000 ms / 512 MB / 1e8 units / 5,000 bp / 500,000 output bp,
# trie 30,000 ms / 512 MB / 1e8 / 1,000 live paths / 2,000 bp): the direct request sends the preset's merged
# strategy, with output.detail "full" (checks.py reads a graphlet body too)
WALK_PRESET = {
    "direction": "both", "exhaustive": False, "labels": {"mode": "constrain"},
    "branching": {"max_label_branches": 2},
    "bounds": {"time_budget_ms": 15000, "max_memory_mb": 512, "max_work_units": 100000000,
               "max_extension_bp": 5000, "max_output_bp": 500000},
}
TRIE_PRESET = {
    "direction": "both", "exhaustive": True, "labels": {"mode": "constrain"},
    "bounds": {"time_budget_ms": 30000, "max_memory_mb": 512, "max_work_units": 100000000,
               "max_live_paths": 1000, "max_extension_bp": 2000},
}
# overrides traverse_fetch can express (DESIGN §4; async's traverse.not_expressible.strategy_keys); the service bounds
# labels.max_seed_labels to its TRAVERSAL_MAX_SEED_LABELS, 10,000 (app/traversal/presets.py SERVICE_BOUNDED: VERIFY-T,
# G14), the host's capabilities state the same maximum
TRAVERSE_OVERRIDABLE = {"bounds", "direction", "exhaustive", "labels.max_labels_per_node", "labels.seed_label_kind",
                        "labels.max_seed_labels"}
SERVER_MAX_SEED_LABELS = 10000
# a label no annotation holds: the builders' unknown label (resolve.unknown_label names it, G10)
ABSENT_LABEL = "NOPE_0000.1"
# The host's substitution caps per pattern kind (/pattern/capabilities mismatches.max), from the capabilities snapshot
# the build reads (set_host_capabilities); the SPEC's maximum until then. pattern_by_design names
# mismatch.above_kind_cap above a kind's cap (G05), and the M builders emit the m = 2 peptides only under a protein cap
# of at least 2 (G04).
HOST_MISMATCH_MAX = {"dna": 3, "iupac": 3, "protein": 3}
# The host's time-budget cap per route from the snapshot (set_host_capabilities): /pattern caps.time_budget_ms,
# /resolve resolve.time_budget.max_time_ms of /traverse/capabilities. A request above it is lowered to it (SPEC-pattern
# §4.5, SPEC-labeled-traversal-core §4.5), so budget_respected states it as cap_ms where a request asks more (A2-03)
HOST_BUDGET_CAP_MS: Dict[str, int] = {}
# The index the snapshot states (/traverse/capabilities index_fp): a claim the records alone cannot certify (a
# graph-wide or column-wide absence) names it as the index it was observed on (areas.empirical)
HOST_INDEX: Dict[str, Optional[str]] = {"index_fp": None}


def _int_or_none(v: Any) -> Optional[int]:
    return v if isinstance(v, int) and not isinstance(v, bool) else None


def set_host_capabilities(pattern_caps: Optional[dict], traverse_caps: Optional[dict] = None) -> None:
    """Take the host's per-kind substitution caps from a /pattern/capabilities document and the routes' time-budget
    caps from it and a /traverse/capabilities document (the build's snapshot)."""
    mx = ((pattern_caps or {}).get("mismatches") or {}).get("max") or {}
    HOST_MISMATCH_MAX.clear()
    HOST_MISMATCH_MAX.update({"dna": 3, "iupac": 3, "protein": 3})
    for kind, v in mx.items():
        if isinstance(v, int) and not isinstance(v, bool):
            HOST_MISMATCH_MAX[kind] = v
    HOST_INDEX["index_fp"] = (traverse_caps or {}).get("index_fp")
    HOST_BUDGET_CAP_MS.clear()
    for route, v in (("POST /pattern", ((pattern_caps or {}).get("caps") or {}).get("time_budget_ms")),
                     ("POST /resolve", (((traverse_caps or {}).get("resolve") or {}).get("time_budget") or {}).get(
                         "max_time_ms"))):
        if _int_or_none(v) is not None:
            HOST_BUDGET_CAP_MS[route] = v

PATTERN_NOT_EXPRESSIBLE = ["max_steps", "max_memory_mb", "min_information_bits", "max_annotation_work",
                           "max_occurrences_per_label", "max_labels_per_anchor", "stop_at_threshold",
                           "allow_unbudgeted_annotation", "graphs", "budget_split"]
PATTERN_MCP_FIELDS = {"patterns", "mode", "scope", "strands", "output", "max_contexts", "time_budget_ms",
                      "long_search", "max_paths", "require_support", "genetic_code", "predicate",
                      "predicate_strands", "predicate_scope", "max_predicate_contexts", "max_predicate_work",
                      "max_labels", "max_anchors", "supported_paths_level", "in_ram", "max_mismatches"}


def deep_merge(base: dict, over: dict) -> dict:
    out = copy.deepcopy(base)
    for k, v in over.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = deep_merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


# the service's trace preset (app/traversal/presets.py: the walk preset with support "trace" and
# branching.on_reconverge "keep", since SPEC-labeled-traversal-core §13 item 9 refuses trace with merging: G15)
TRACE_PRESET = deep_merge(WALK_PRESET, {"support": "trace", "branching": {"on_reconverge": "keep"}})


def canonical_json(obj: Any) -> str:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(canonical_json(obj).encode("utf-8")).hexdigest()


class Ctx:
    """The records with sequences and the seeded random streams, one per (area, purpose)."""

    def __init__(self, records: List[Record], seed: int):
        self.records = records
        self.by_key = {r.key: r for r in records}
        self.seed = seed
        self.positives = [r for r in records if r.usable_for_positives and r.key != "athc"]
        self._used: Dict[str, set] = {}
        self.relations: List[dict] = []   # the manifest's cross-scenario comparisons (relation())
        self.skipped: List[dict] = []     # planned scenarios a generator could not build, with the reason

    def skip(self, area: str, what: str, reason: str) -> None:
        """State a planned scenario that could not be built (the manifest lists them; nothing drops silently)."""
        self.skipped.append({"area": area, "scenario": what, "reason": reason})

    def rng(self, *parts: Any) -> random.Random:
        h = hashlib.sha256(("%d:" % self.seed + ":".join(str(p) for p in parts)).encode()).hexdigest()
        return random.Random(int(h[:16], 16))

    def rec(self, key: str) -> Record:
        return self.by_key[key]

    # ------------------------------------------------------------------------------------------ windows
    def window(self, rng: random.Random, rec: Record, L: int, margin: int = 40, unique: bool = True,
               max_dust: float = 2.0, where: Optional[Tuple[int, int]] = None, tries: int = 2000) -> dict:
        """A window of L bases of |rec| (ACGT only, low DUST score, unique on both strands of the fetched span
        when |unique|), at least |margin| bases from the span's ends; |where| restricts the start to a record
        coordinate range. Returns the truth fields."""
        seq = rec.seq
        lo, hi = margin, len(seq) - L - margin
        if where:
            lo = max(lo, where[0] - rec.offset)
            hi = min(hi, where[1] - rec.offset)
        if hi < lo:
            raise ValueError("no room for a %d bp window in %s" % (L, rec.key))
        used = self._used.setdefault(rec.key, set())
        for _ in range(tries):
            s = rng.randint(lo, hi)
            w = seq[s:s + L]
            if not st.is_acgt(w) or (s // 50) in used:
                continue
            if st.low_complexity_dust_score(w[:64]) > max_dust:
                continue
            if unique:
                if seq.count(w) != 1 or (st.revcomp(w) != w and st.revcomp(w) in seq):
                    continue
            used.add(s // 50)
            return self.truth_at(rec, s, L)
        raise ValueError("no window of %d bp found in %s" % (L, rec.key))

    def truth_at(self, rec: Record, s: int, L: int, flank: int = K - 1, strand: str = "+") -> dict:
        """Truth of the bases rec.seq[s:s+L] (s an index into the fetched span): record coordinates (0-based start),
        the window as the record holds it on its forward strand, and a flank to recompute the source contexts (cut at
        a base outside ACGT: acgt_flank)."""
        seq = rec.seq
        a, b, cut = acgt_flank(rec, s, L, flank)
        out = {
            "record": rec.accession, "genome": rec.key, "taxid": rec.column, "seq_id": rec.seq_id,
            "record_length": rec.length, "start": rec.offset + s, "length": L, "strand": strand,
            "window": seq[s:s + L], "flank": seq[a:b], "flank_start": rec.offset + a,
            "fetched_span": [rec.offset, rec.offset + len(seq)],
        }
        if cut:
            out["flank_cut"] = cut
        return out

    def rc_negative_window(self, rng: random.Random, rec: Record, L: int, tries: int = 50, **kw: Any) -> dict:
        """A window as window() picks it whose reverse complement shares no k-mer with the record's forward strand
        (the whole record), so that the record's header supports none of the reverse query's k-mers on this basic
        graph: whole-window uniqueness on both strands does not imply it (a window overlapping an inverted repeat keeps
        reverse k-mers in the record, about 1 % of the windows window() accepts). A window drawn after a rejected one
        comes from the same stream. Raises when none is found."""
        for _ in range(tries):
            t = self.window(rng, rec, L, **kw)
            if not reverse_kmers_in_record(rec, t["window"]):
                return t
        raise ValueError("no %d bp window of %s without a reverse k-mer in its record" % (L, rec.key))

    def flank_for(self, rec: Record, t: dict, extra: int) -> dict:
        """The same truth with a flank of |extra| bases each side (traversal walks), cut at a base outside ACGT
        (acgt_flank)."""
        s = t["start"] - rec.offset
        out = dict(t)
        a, b, cut = acgt_flank(rec, s, t["length"], extra)
        out["flank"] = rec.seq[a:b]
        out["flank_start"] = rec.offset + a
        out.pop("flank_cut", None)
        if cut:
            out["flank_cut"] = cut
        return out


def acgt_flank(rec: Record, s: int, L: int, extra: int) -> Tuple[int, int, Dict[str, dict]]:
    """The flank [a, b) of the window rec.seq[s:s+L] (indices into the fetched span), |extra| bases each side within
    the span, cut before the first base outside ACGT on each side: no k-mer of the graph holds such a base (the
    alphabet is $ACGT, SPEC-pattern-search §8.2), so a walk cannot pass it and no source context holds it. The cut is
    stated per side as {at (the record coordinate of that base), base}: a structural end of the truth like the
    record's end (T-walk-0018's left flank held an IUPAC S 4,898 bp from its seed, triage-v4 V03)."""
    seq = rec.seq
    a = max(0, s - extra)
    b = min(len(seq), s + L + extra)
    cut: Dict[str, dict] = {}
    for x in range(s - 1, a - 1, -1):
        if seq[x] not in "ACGT":
            cut["left"] = {"at": rec.offset + x, "base": seq[x]}
            a = x + 1
            break
    for x in range(s + L, b):
        if seq[x] not in "ACGT":
            cut["right"] = {"at": rec.offset + x, "base": seq[x]}
            b = x
            break
    return a, b, cut


def reverse_kmers_in_record(rec: Record, window: str, k: int = K) -> List[int]:
    """The starts of the k-mers of revcomp(window) that the record's forward strand holds. It needs the whole record:
    a fetched region cannot certify the rest of the record."""
    if rec.region:
        raise ValueError("%s: only a region is fetched; its reverse k-mers cannot be certified absent" % rec.key)
    q = st.revcomp(window)
    return [i for i in range(len(q) - k + 1) if q[i:i + k] in rec.seq]


# ---------------------------------------------------------------------------------------------- walk expectations
def walk_expectations(seq: str, plans: List[Tuple[int, int, List[str], int]], k: int = K) -> List[Dict[str, dict]]:
    """For walks constrained to the header of the record |seq| (whole, both strands considered), per plan
    (s0 = seed start in seq, L = seed length, sides, E = max_extension_bp) and side: {"min_bp", "first",
    "branch_by_bp"}.

    min_bp: how far the walk must follow the record before any stated end is legitimate: the first step whose
    (k-1)-mer context has a second successor (right; a second predecessor, left) anywhere in the record on either
    strand (the header may then label a branch), the first edge (k+1-mer) or seed k-mer the walk meets again on
    either strand (edge_reuse / rejoined_seed), the record's end, or E. branch_by_bp: the first step whose
    context has a second successor on the forward strand (the record's own copies, which the header certainly
    labels: the walk must state a split or branch event by then), None when there is none within E.
    A base outside ACGT ends the record's walk where the next step would take it ("first": "non_acgt"; no k-mer of
    the graph holds it), and a k-mer holding one is no successor anywhere.
    One pass per strand over the record for all plans (the verification's h1c.py method)."""
    need: Dict[Tuple[str, str], List[set]] = {}
    for s0, L, sides, E in plans:
        for side in sides:
            for j in range(E):
                if side == "right":
                    p = s0 + L + j
                    if p >= len(seq):
                        break
                    need.setdefault(("R", seq[p - k + 1:p]), [set(), set()])
                else:
                    p = s0 - 1 - j
                    if p < 0:
                        break
                    need.setdefault(("L", seq[p + 1:p + k]), [set(), set()])
    ctx_r = {c: v for (d, c), v in need.items() if d == "R"}
    ctx_l = {c: v for (d, c), v in need.items() if d == "L"}
    for strand, text in ((0, seq), (1, st.revcomp(seq))):
        for i in range(len(text) - k + 1):
            km = text[i:i + k]
            if km[0] not in "ACGT" or km[-1] not in "ACGT":
                continue   # no k-mer of the graph (a context looked up is ACGT, so only its added base can be another)
            v = ctx_r.get(km[:-1])
            if v is not None:
                v[0].add(km[-1])
                if strand == 0:
                    v[1].add(km[-1])
            v = ctx_l.get(km[1:])
            if v is not None:
                v[0].add(km[0])
                if strand == 0:
                    v[1].add(km[0])
    out = []
    for s0, L, sides, E in plans:
        res: Dict[str, dict] = {}
        seed = seq[s0:s0 + L]
        seen_edges = set()
        for x in range(L - k):
            e = seed[x:x + k + 1]
            seen_edges |= {e, st.revcomp(e)}
        seed_kmers = set()
        for x in range(L - k + 1):
            km = seed[x:x + k]
            seed_kmers |= {km, st.revcomp(km)}
        for side in sides:
            edges = set(seen_edges)
            first, at, branch_by = "max_extension", E, None
            for j in range(E):
                if side == "right":
                    p = s0 + L + j
                    if p >= len(seq):
                        first, at = "record_end", j
                        break
                    if seq[p] not in "ACGT":
                        first, at = "non_acgt", j
                        break
                    succ = ctx_r[seq[p - k + 1:p]]
                    edge, kmer = seq[p - k:p + 1], seq[p - k + 1:p + 1]
                else:
                    p = s0 - 1 - j
                    if p < 0:
                        first, at = "record_end", j
                        break
                    if seq[p] not in "ACGT":
                        first, at = "non_acgt", j
                        break
                    succ = ctx_l[seq[p + 1:p + k]]
                    edge, kmer = seq[p:p + k + 1], seq[p:p + k]
                if branch_by is None and len(succ[1]) > 1:
                    branch_by = j
                if len(succ[0]) > 1:
                    first, at = "ambiguous", j
                    break
                if edge in edges or kmer in seed_kmers:
                    first, at = "reentry", j
                    break
                edges |= {edge, st.revcomp(edge)}
            if branch_by is None:
                # the forward-strand scan may reach past the first both-strand ambiguity
                for j in range(at, E):
                    if side == "right":
                        p = s0 + L + j
                        if p >= len(seq) or seq[p] not in "ACGT":
                            break
                        succ = ctx_r.get(seq[p - k + 1:p])
                    else:
                        p = s0 - 1 - j
                        if p < 0 or seq[p] not in "ACGT":
                            break
                        succ = ctx_l.get(seq[p + 1:p + k])
                    if succ is not None and len(succ[1]) > 1:
                        branch_by = j
                        break
            res[side] = {"min_bp": at, "first": first, "branch_by_bp": branch_by}
        out.append(res)
    return out


# ---------------------------------------------------------------------------------------------- requests
def pattern_request(patterns: List[dict], mode: str = "all_or_count", max_contexts: int = PATTERN_MAX_CONTEXTS,
                    time_budget_ms: int = PATTERN_TIME_BUDGET_MS, **fields: Any) -> dict:
    """A /pattern body that names mode, max_contexts and time_budget_ms (async finding 1) and the defaults both arms
    must agree on wherever they act (the contract: one question on every arm): scope, strands and output.labels
    always; output.occurrences with labels "all" or "predicate_only"; predicate_strands with a predicate. Where
    they do not act they stay unnamed: occurrences with labels "none" is a host 400 when true and a service
    refusal for any value (async: pattern.occurrences_without_labels), predicate_strands without a predicate a
    service refusal (pattern.predicate_options_without_predicate) and a no-op on the host."""
    req = {"patterns": patterns, "mode": mode, "max_contexts": max_contexts, "time_budget_ms": time_budget_ms}
    for k, v in fields.items():
        if v is not None:
            req[k] = v
    req.setdefault("scope", PATTERN_DEFAULTS["scope"])
    req.setdefault("strands", PATTERN_DEFAULTS["strands"])
    out = dict(req.get("output") or {})
    out.setdefault("labels", PATTERN_DEFAULTS["labels"])
    if out["labels"] in ("all", "predicate_only"):
        out.setdefault("occurrences", PATTERN_DEFAULTS["occurrences"])
    req["output"] = out
    if req.get("predicate") is not None:
        req.setdefault("predicate_strands", PATTERN_DEFAULTS["predicate_strands"])
    if req.get("output", {}).get("labels") in ("all", "predicate_only") and mode == "partial" \
            and "max_labels" not in req:
        req["max_labels"] = 1000   # the host's default, named: the MCP's default is 100
    if req.get("predicate") is not None and "max_predicate_contexts" not in req:
        req["max_predicate_contexts"] = PATTERN_MAX_PREDICATE_CONTEXTS   # the MCP's default, named on both arms
    return req


def resolve_request(sequence: str, labels: Optional[List[str]] = None, discover: Optional[dict] = None,
                    time_budget_ms: Optional[float] = RESOLVE_TIME_BUDGET_MS, max_query_bp: Optional[int] = None,
                    **fields: Any) -> dict:
    req: Dict[str, Any] = {"sequence": sequence}
    if labels is not None:
        req["labels"] = labels
    if discover is not None:
        req["discover"] = discover
    bounds: Dict[str, Any] = {}
    if max_query_bp is None:
        max_query_bp = DISCOVER_MAX_QUERY_BP if discover is not None else RESOLVE_MAX_QUERY_BP
    if max_query_bp and len(sequence) <= max_query_bp:
        bounds["max_query_bp"] = max_query_bp
    if time_budget_ms is not None and discover is None:
        bounds["time_budget_ms"] = time_budget_ms
    if bounds:
        req["bounds"] = bounds
    req.update(fields)
    return req


def traverse_request(seeds: List[dict], overrides: Optional[dict] = None, preset: str = "walk",
                     detail: str = "full", extra: Optional[dict] = None) -> dict:
    """A /traverse body: the MCP preset's merged strategy (walk, trie or trace), then the overrides, then |extra|
    (fields the MCP cannot express), with output.detail."""
    base = {"walk": WALK_PRESET, "trie": TRIE_PRESET, "trace": TRACE_PRESET}[preset]
    strategy = deep_merge(base, overrides or {})
    strategy = deep_merge(strategy, {"output": {"detail": detail, "timing": True}})
    if extra:
        strategy = deep_merge(strategy, extra)
    return {"seeds": seeds, "strategy": strategy}


# ---------------------------------------------------------------------------------------------- MCP mapping
_STRANDS_MCP = {"both": "both", "forward": "plus", "reverse": "minus"}


def _flat_keys(d: dict, prefix: str = "") -> List[str]:
    out = []
    for k, v in d.items():
        p = prefix + k
        if isinstance(v, dict):
            out += _flat_keys(v, p + ".")
        else:
            out.append(p)
    return out


def pattern_by_design(req: Any) -> List[str]:
    """The async table's case ids a /pattern request matches (the rules its scenario_match states)."""
    ids: List[str] = []
    if not isinstance(req, dict):
        return ids
    pats = req.get("patterns")
    if isinstance(pats, list):
        if len(pats) == 0:
            ids.append("pattern.empty_pattern_list")
        elif 16 < len(pats) <= 64:
            ids.append("pattern.too_many_patterns_chunked")
        elif len(pats) > 64:
            ids.append("pattern.over_64_patterns")
        kinds = set()
        idlist = []
        for p in pats:
            if not isinstance(p, dict):
                continue
            named = [k for k in ("dna", "iupac", "protein") if k in p]
            if len(named) != 1:
                ids.append("pattern.two_kinds_in_one_pattern")
                continue
            kind = named[0]
            kinds.add("protein" if kind == "protein" else "nucleotide")
            s = p[kind]
            if "id" in p:
                idlist.append(p["id"])
            if not isinstance(s, str):
                continue
            if s == "":
                ids.append("pattern.empty_pattern_string")
            if len(s) > 1000:
                ids.append("pattern.length_above_1000")
            if s != s.upper():
                ids.append("pattern.lowercase")
            up = s.upper()
            if kind in ("dna", "iupac") and "U" in up:
                ids.append("pattern.u_for_t")
            elif kind == "dna" and any(c not in "ACGT" for c in up) and up:
                ids.append("pattern.bad_alphabet_dna")
            elif kind == "iupac" and any(c not in st.IUPAC for c in up):
                ids.append("pattern.bad_alphabet_iupac")
            elif kind == "protein" and any(c not in st.AMINO + "XBZJ*" for c in up):
                ids.append("pattern.protein_alphabet")
            L = len(s) * (3 if kind == "protein" else 1)
            if L > K:
                ls = req.get("long_search")
                if req.get("predicate") is not None and ls in (None, "anchors"):
                    ids.append("pattern.predicate_long_pattern_without_supported_paths")
                elif ls is None and req.get("predicate") is None:
                    ids.append("pattern.long_without_long_search")
                if req.get("max_mismatches") not in (None, 0):
                    ids.append("mismatch.long_pattern_with_long_search" if ls == "paths"
                               else "mismatch.long_pattern_without_long_search")
        if len(kinds) > 1:
            ids.append("pattern.mixed_nucleotide_protein")
        if len(idlist) != len(set(map(str, idlist))):
            ids.append("pattern.duplicate_ids")
        if req.get("genetic_code") is not None and "protein" not in kinds:
            ids.append("pattern.genetic_code_without_protein")
    if req.get("strands") in ("forward", "reverse"):
        ids.append("pattern.strands_vocabulary")
    tb = req.get("time_budget_ms")
    if isinstance(tb, (int, float)) and not isinstance(tb, bool):
        if tb > 600000:
            ids.append("pattern.time_budget_above_host_cap")
        elif tb < 500:
            ids.append("pattern.time_budget_below_minimum")
    mc = req.get("max_contexts")
    if isinstance(mc, int) and mc > 10000:
        ids.append("pattern.max_contexts_above_10000")
    # a cap of 0, which the host answers (SPEC-pattern §4.1: integers >= 0) and the service's REST model refuses
    # (caps >= 1; async's case of G02)
    if any(req.get(f) == 0 and not isinstance(req.get(f), bool) for f in ("max_contexts", "max_anchors", "max_paths")):
        ids.append("pattern.caps_zero")
    labels = (req.get("output") or {}).get("labels", "none")
    mode = req.get("mode", "all_or_count")
    if "max_labels" in req and (labels == "none" or mode != "partial"):
        ids.append("pattern.max_labels_wrong_mode_or_projection")
    if "occurrences" in (req.get("output") or {}) and labels == "none":
        ids.append("pattern.occurrences_without_labels")
    ls = req.get("long_search")
    if ls == "anchors":
        ids.append("pattern.long_search_anchors")
    if ls is None and ("max_paths" in req or "max_anchors" in req):
        ids.append("pattern.path_caps_without_long_search")
    if "supported_paths_level" in req and ls != "supported_paths":
        ids.append("pattern.supported_paths_level_without_supported_paths")
    gc = req.get("genetic_code")
    if isinstance(gc, int) and gc not in st.NCBIEAA:
        ids.append("pattern.genetic_code_unknown")
    pred = req.get("predicate")
    if pred is not None:
        if ls == "paths":
            ids.append("pattern.predicate_with_paths")
        if "max_predicate_contexts" in req and isinstance(req["max_predicate_contexts"], int) \
                and req["max_predicate_contexts"] > 10000:
            ids.append("pattern.max_predicate_contexts_above_10000")
        if _predicate_has_number(pred):
            ids.append("pattern.predicate_name_not_string")
        if _predicate_form_fault(pred):
            ids.append("pattern.predicate_form_other")
    else:
        if labels == "predicate_only":
            ids.append("pattern.predicate_only_without_predicate")
        # split by value in async's table (789fb5f0, G08): "motif" without a predicate is a host refusal (the old
        # id, what suite_v3's Q-motif-0004 sends), predicate_strands or "context" alone a host no-op
        if req.get("predicate_scope") == "motif":
            ids.append("pattern.predicate_options_without_predicate")
        elif "predicate_strands" in req or "predicate_scope" in req:
            ids.append("pattern.predicate_options_noop_without_predicate")
        if "max_predicate_contexts" in req or "max_predicate_work" in req:
            ids.append("pattern.predicate_caps_without_predicate")
    if "max_mismatches" in req:
        m = req["max_mismatches"]
        if m is None:
            ids.append("mismatch.value_null")
        elif isinstance(m, bool):
            ids.append("mismatch.value_true")
        elif isinstance(m, str):
            ids.append("mismatch.value_string_1")
        elif isinstance(m, float) and not m.is_integer():
            ids.append("mismatch.value_1_5")
        elif m == 0:
            ids.append("mismatch.zero")
        elif m < 0:
            ids.append("mismatch.value_minus_1")
        elif m > 3:
            ids.append("mismatch.value_4")
        elif m > 0:
            # above a kind's cap of the host (staging: protein 1), for a pattern of at most k bases (a longer one is
            # the long_pattern case): the service refuses the whole submit, the host answers that pattern's slot
            if any(m > HOST_MISMATCH_MAX.get(kind, 3) for kind, n in _pattern_kinds(pats) if n <= K):
                ids.append("mismatch.above_kind_cap")
            if labels in ("all", "predicate_only"):
                ids.append("mismatch.with_labels")
            if pred is not None and labels == "none":
                ids.append("mismatch.with_predicate")
            if ls == "supported_paths":
                ids.append("mismatch.with_supported_paths")
    for f in PATTERN_NOT_EXPRESSIBLE:
        if f in req:
            ids.append("pattern.not_expressible." + f)
    if "paths" in (req.get("output") or {}):
        ids.append("pattern.not_expressible.output_paths")
    return ids


def _pattern_kinds(pats: Any) -> List[Tuple[str, int]]:
    """(kind, length in bases) of each pattern of a request (a residue counts three bases)."""
    out = []
    for p in pats if isinstance(pats, list) else []:
        if not isinstance(p, dict):
            continue
        named = [k for k in ("dna", "iupac", "protein") if k in p]
        if len(named) == 1 and isinstance(p[named[0]], str):
            out.append((named[0], len(p[named[0]]) * (3 if named[0] == "protein" else 1)))
    return out


_PREDICATE_OPERATORS = ("any", "all", "none", "at_least", "and", "or", "not")


def _predicate_form_fault(p: Any, depth: int = 0) -> bool:
    """Whether a predicate breaks a form rule of SPEC-pattern-search §19.3 other than a non-string name (that is
    pattern.predicate_name_not_string): an object of not exactly one operator, an unknown operator, an empty name
    list or a name twice in one, an empty name, at_least with other keys or n outside 1..len, an empty and/or list,
    nesting deeper than 64 (async's case pattern.predicate_form_other)."""
    if depth > 64 or not isinstance(p, dict) or len(p) != 1:
        return True
    op, arg = next(iter(p.items()))
    if op not in _PREDICATE_OPERATORS:
        return True

    def names_fault(lst: Any) -> bool:
        if not isinstance(lst, list) or not lst:
            return True
        strs = [x for x in lst if isinstance(x, str)]
        return any(x == "" for x in strs) or len(strs) != len(set(strs))

    if op in ("any", "all", "none"):
        return names_fault(arg)
    if op == "at_least":
        if not isinstance(arg, dict) or set(arg) != {"n", "labels"} or names_fault(arg["labels"]):
            return True
        n = arg["n"]
        return not (isinstance(n, int) and not isinstance(n, bool) and 1 <= n <= len(arg["labels"]))
    if op in ("and", "or"):
        return not isinstance(arg, list) or not arg or any(_predicate_form_fault(x, depth + 1) for x in arg)
    return _predicate_form_fault(arg, depth + 1)


def _predicate_has_number(p: Any) -> bool:
    if isinstance(p, dict):
        return any(_predicate_has_number(v) for v in p.values())
    if isinstance(p, list):
        return any(isinstance(x, (int, float)) and not isinstance(x, bool) or _predicate_has_number(x) for x in p)
    return False


def mcp_pattern(req: Any) -> dict:
    """pattern_search arguments for a /pattern body (DESIGN §4), or not_implemented naming the field."""
    by_design = pattern_by_design(req)
    if not isinstance(req, dict):
        return {"not_implemented": "non-object body", "by_design": ["capabilities.not_expressible.protocol"]}
    for f in PATTERN_NOT_EXPRESSIBLE:
        if f in req:
            return {"not_implemented": f, "by_design": by_design}
    unknown = [k for k in req if k not in PATTERN_MCP_FIELDS]
    if unknown:
        return {"not_implemented": unknown[0], "by_design": by_design + ["pattern.unknown_field"]}
    out_block = req.get("output") or {}
    if "paths" in out_block:
        return {"not_implemented": "output.paths", "by_design": by_design}
    # the MCP coerces integers ("1", true, 1.0 -> 1) and drops what it has no field for (async finding 2): a
    # host refusal of such a value cannot be asked through it
    m = req.get("max_mismatches", 0)
    if "max_mismatches" in req and (m is None or isinstance(m, (bool, str)) or (isinstance(m, float))):
        return {"not_implemented": "max_mismatches (value %s: the MCP coerces or refuses it)" % json.dumps(m),
                "by_design": by_design}
    args: Dict[str, Any] = {"databases": [DB], "patterns": copy.deepcopy(req.get("patterns", [])),
                            "mode": req.get("mode", "all_or_count"),
                            "scope": req.get("scope", "any_offset"),
                            "strands": _STRANDS_MCP.get(req.get("strands", "both"), req.get("strands")),
                            "labels": out_block.get("labels", "none")}
    if "occurrences" in out_block:
        args["occurrences"] = out_block["occurrences"]
    for f in ("max_contexts", "time_budget_ms", "long_search", "max_paths", "require_support", "genetic_code",
              "predicate", "predicate_strands", "predicate_scope", "max_predicate_contexts",
              "max_predicate_work", "max_labels", "max_anchors", "supported_paths_level", "max_mismatches"):
        if f in req:
            args[f] = copy.deepcopy(req[f])
    if "in_ram" in req:
        # the service's documented mapping: high_throughput sends in_ram (async: pattern.in_ram_by_size)
        args["high_throughput"] = req["in_ram"]
        by_design = sorted(set(by_design) | {"pattern.in_ram_by_size"})
    block = {"tool": "pattern_search", "arguments": args}
    if req.get("max_mismatches", 0):
        block["requires_service"] = "I12"
    if by_design:
        block["by_design"] = by_design
    return block


def resolve_by_design(req: Any) -> List[str]:
    ids = ["resolve.job_path"]
    if not isinstance(req, dict):
        return ids
    seq = req.get("sequence", "")
    labels = req.get("labels")
    disc = req.get("discover")
    if isinstance(seq, str):
        if seq.strip() == "":
            ids.append("resolve.empty_sequence")
        else:
            if len(seq) < K:
                ids.append("resolve.sequence_too_short")
            if labels is not None and len(seq) > 2000:
                ids.append("resolve.sequence_too_long")
            if disc is not None and 500 < len(seq) <= 2000:
                ids.append("resolve.discovery_too_long")
            up = seq.upper()
            if any(c not in "ACGT" for c in up):
                ids.append("resolve.n_or_iupac")
            if seq != up:
                ids.append("resolve.lowercase_whitespace_fasta")
    if labels is not None and disc is not None:
        ids.append("resolve.discover_with_labels")
    if labels is None and disc is None:
        ids.append("resolve.no_labels_no_discover")
    if isinstance(labels, list):
        if ABSENT_LABEL in labels:
            ids.append("resolve.unknown_label")
        if len(labels) != len(set(labels)):
            ids.append("resolve.duplicate_labels")
        if len(set(labels)) > 1000:
            ids.append("resolve.too_many_labels")
        if any(isinstance(l, str) and l.strip() == "" for l in labels):
            ids.append("resolve.blank_label")
    # split by value in async's table (789fb5f0, G11): the old ids keep the value suite_v3 sends (max_labels 0, a
    # budget of 250 ms, both host refusals); above 100 labels and 251..499 ms the host answers
    ml = disc.get("max_labels") if isinstance(disc, dict) else None
    if isinstance(ml, int) and not isinstance(ml, bool):
        if ml < 1:
            ids.append("resolve.discover_limit_above_100")
        elif ml > 100:
            ids.append("resolve.discover_limit_above_100_host_answers")
    tb = (req.get("bounds") or {}).get("time_budget_ms")
    if isinstance(tb, (int, float)) and not isinstance(tb, bool):
        if tb > HOST_BUDGET_CAP_MS.get("POST /resolve", 60000):
            ids.append("resolve.time_budget_above_host_max")
        elif tb <= 250:
            ids.append("resolve.time_budget_below_minimum")
        elif tb < 500:
            ids.append("resolve.time_budget_above_reserve_below_minimum")
    elif labels is not None and disc is None:
        ids.append("resolve.time_budget_unnamed")
    mq = (req.get("bounds") or {}).get("max_query_bp")
    if mq != (DISCOVER_MAX_QUERY_BP if disc is not None else RESOLVE_MAX_QUERY_BP):
        ids.append("resolve.not_expressible.bounds_max_query_bp")
    for f in ("support", "select", "graph"):
        if f in req:
            ids.append("resolve.not_expressible." + f)
    return ids


def mcp_resolve(req: Any) -> dict:
    by_design = resolve_by_design(req)
    if not isinstance(req, dict):
        return {"not_implemented": "non-object body", "by_design": ["capabilities.not_expressible.protocol"]}
    for f in ("support", "select", "graph", "graphs", "min_block_kmers", "run_format"):
        if f in req:
            return {"not_implemented": f, "by_design": by_design}
    known = {"sequence", "labels", "discover", "bounds"}
    unknown = [k for k in req if k not in known]
    if unknown:
        return {"not_implemented": unknown[0], "by_design": by_design}
    args: Dict[str, Any] = {"database": DB}
    if "sequence" in req:
        args["sequence"] = req["sequence"]
    if "labels" in req:
        args["labels"] = req["labels"]
    if "discover" in req:
        args["discover"] = True
        d = req["discover"] if isinstance(req["discover"], dict) else {}
        if "max_labels" in d:
            args["limit"] = d["max_labels"]
        if "kind" in d:
            args["kind"] = d["kind"]
    b = req.get("bounds") or {}
    if "time_budget_ms" in b:
        args["time_budget_ms"] = b["time_budget_ms"]
    block = {"tool": "traverse_resolve", "arguments": args}
    mq = b.get("max_query_bp")
    if mq is not None and mq != (DISCOVER_MAX_QUERY_BP if "discover" in req else RESOLVE_MAX_QUERY_BP):
        return {"not_implemented": "bounds.max_query_bp", "by_design": by_design}
    block["by_design"] = by_design
    return block


def traverse_by_design(req: Any) -> List[str]:
    ids: List[str] = []
    if not isinstance(req, dict):
        return ["traverse.not_expressible.non_object_body"]
    seeds = req.get("seeds") or []
    if isinstance(req.get("seeds"), list) and not req["seeds"]:
        ids.append("traverse.empty_seeds")   # traverse_fetch refuses it locally (async's case of G19)
    if isinstance(seeds, list):
        seqs = []
        for s in seeds:
            if not isinstance(s, dict):
                continue
            q = s.get("sequence", "")
            seqs.append((q, tuple(s.get("labels") or [])))
            if isinstance(q, str):
                if len(q) < K:
                    ids.append("traverse.seed_shorter_than_k")
                if any(c not in "ACGT" for c in q.upper()):
                    ids.append("traverse.seed_with_n_or_iupac")
                if len(q) > 55000:
                    ids.append("traverse.seed_too_long")
            if "labels" in s and s["labels"] == []:
                ids.append("traverse.empty_seed_labels")
        if 4 < len(seeds) <= 10000:
            ids.append("traverse.seeds_above_host_max_seeds")
        if len(seqs) != len(set(seqs)):
            ids.append("traverse.duplicate_seeds")
    strategy = req.get("strategy") or {}
    preset = TRIE_PRESET if strategy.get("exhaustive") else WALK_PRESET
    for key in _flat_keys(strategy):
        top = key.split(".")[0]
        if key.startswith("output.") or key == "schema_version":
            if key not in ("output.detail", "output.timing"):
                ids.append("traverse.not_expressible.strategy_keys")
            continue
        if top not in ("direction", "support", "exhaustive", "labels", "branching", "bounds", "frontier",
                       "annotation", "output", "schema_version"):
            ids.append("traverse.unknown_strategy_key")
    b = strategy.get("bounds") or {}
    for knob in b:
        if knob not in ("max_extension_bp", "min_live_labels", "max_steps", "max_live_paths", "max_paths",
                        "max_output_bp", "time_budget_ms", "max_memory_mb", "max_work_units"):
            ids.append("traverse.unknown_bounds_key")
    row = preset["bounds"]
    for knob in ("time_budget_ms", "max_memory_mb", "max_work_units"):
        if isinstance(b.get(knob), (int, float)) and b[knob] > row.get(knob, float("inf")):
            ids.append("traverse.bounds_budget_above_preset_row")
    maxima = {"max_extension_bp": 10000, "max_steps": 1000000, "max_paths": 10000, "max_output_bp": 2000000,
              "max_live_paths": 1000}
    if any(isinstance(b.get(k), int) and b[k] > v for k, v in maxima.items()):
        ids.append("traverse.bounds_knob_above_service_maximum")
    if strategy.get("support") == "trace":
        ids.append("traverse.trace_preset")
    for f in ("attempt_id", "budget_id", "locus_id", "not_after_ms", "expect_server_instance"):
        if f in req:
            ids.append("traverse.not_expressible.attempt_fields")
            break
    for f in ("graph", "release"):
        if f in req:
            ids.append("traverse.not_expressible." + f)
    return ids


def mcp_traverse(req: Any) -> dict:
    """traverse_fetch arguments: the preset (trie when exhaustive), the overrides it can express; a strategy key
    outside them that differs from the preset's value is not implemented yet via MCP."""
    by_design = traverse_by_design(req)
    if not isinstance(req, dict):
        return {"not_implemented": "non-object body", "by_design": by_design}
    for f in ("attempt_id", "budget_id", "locus_id", "not_after_ms", "expect_server_instance", "graph",
              "graphs", "release", "in_ram"):
        if f in req:
            return {"not_implemented": f, "by_design": by_design}
    strategy = req.get("strategy") or {}
    preset_name = "trie" if strategy.get("exhaustive") else "walk"
    preset = TRIE_PRESET if preset_name == "trie" else WALK_PRESET
    if strategy.get("support") == "trace":
        preset_name, preset = "trace", TRACE_PRESET
    overrides: Dict[str, Any] = {}
    for key in _flat_keys(strategy):
        parts = key.split(".")
        val = strategy
        for p in parts:
            val = val[p]
        pval: Any = preset
        for p in parts:
            pval = pval.get(p) if isinstance(pval, dict) else None
        if key in ("output.detail", "output.timing", "schema_version"):
            continue   # the service exports the graphlet; the detail is the arm's own
        if key == "support" and val == "trace":
            continue   # preset trace (async: traverse.trace_preset)
        if pval == val and key != "bounds.max_extension_bp":
            continue   # the preset's own value
        if parts[0] in ("bounds", "direction", "exhaustive") or key in TRAVERSE_OVERRIDABLE:
            cur = overrides
            for p in parts[:-1]:
                cur = cur.setdefault(p, {})
            cur[parts[-1]] = copy.deepcopy(val)
        else:
            extra_ids = ["traverse.not_expressible.strategy_keys"]
            if key == "labels.mode" and val == "annotate":
                extra_ids.append("traverse.not_expressible.advanced_strategy")
            return {"not_implemented": "strategy." + key,
                    "by_design": sorted(set(by_design + extra_ids))}
    seeds = []
    for s in req.get("seeds") or []:
        if isinstance(s, dict):
            seeds.append({k: copy.deepcopy(v) for k, v in s.items() if k in ("sequence", "labels", "seed_id")})
        else:
            seeds.append(s)
    args = {"database": DB, "seeds": seeds, "preset": preset_name, "force_fresh": True}
    if overrides:
        args["overrides"] = overrides
    block = {"tool": "traverse_fetch", "arguments": args}
    if by_design:
        block["by_design"] = by_design
    return block


# ---------------------------------------------------------------------------------------------- the time budget
# The rule budget_respected applies per route (stated in the check's args): no fixed allowance, since the
# capabilities state no bound on one uninterruptible piece (deadline_check.max_uninterruptible_ms null)
BUDGET_RULE = {
    "POST /traverse": "SPEC-labeled-traversal-core.md §6.8 and walker.cpp: the deadline is the budget on the seed's "
                      "clock, never read at depth 0, the depth-0 reads unpaced; a time stop not before the budget; "
                      "a budget that passed in the seed phase stops the walk at the first depth-1 checkpoint (no path "
                      "beyond 1 bp); a stop at a checked depth (at_bp >= 2) within one read and one head piece after "
                      "max(budget, seed phase), at depth 1 within the depth-0 levels' pieces; no time stop only "
                      "within the pieces before a reading that would stop a head below the radius",
    "POST /pattern": "SPEC-pattern-search.md §7.6: an answer that cannot be written by time_budget_ms is a 503 "
                     "deadline, so a 200's timing.elapsed_ms (to the end of the search) is within the budget; an "
                     "allowed 503 deadline long past the budget (client wall time) is a finding",
    "POST /resolve": "SPEC-labeled-traversal-core.md §4.5: past the budget the answer is a 503, never a partial "
                     "answer, so a 200's timing.elapsed_ms is within the budget; an allowed 503 deadline long past "
                     "the budget (client wall time) is a finding",
}


def time_budget_of(route: str, req: Any) -> Optional[Any]:
    """The time budget a request names (None without one)."""
    if not isinstance(req, dict):
        return None
    if route == "POST /pattern":
        return req.get("time_budget_ms")
    if route == "POST /resolve":
        return (req.get("bounds") or {}).get("time_budget_ms")
    if route == "POST /traverse":
        return ((req.get("strategy") or {}).get("bounds") or {}).get("time_budget_ms")
    return None


def budget_expectation(route: str, req: Any, expect: List[dict]) -> Optional[dict]:
    """budget_respected for a scenario that names a time budget and allows a 200 answer; with cap_ms, the route's
    cap of the snapshot (HOST_BUDGET_CAP_MS), where the request asks more: the effective budget an allowed 503 is
    judged against is min(requested, cap) (A2-03; below the cap the requested budget is the effective one)."""
    asked = time_budget_of(route, req)
    if asked is None or route not in BUDGET_RULE:
        return None
    allowed = [200]
    for e in expect:
        if e.get("check") == "status":
            w = e.get("args", {}).get("http", 200)
            allowed = w if isinstance(w, list) else [w]
    if 200 not in allowed:
        return None
    args = {"rule": BUDGET_RULE[route]}
    cap = HOST_BUDGET_CAP_MS.get(route)
    if cap is not None and isinstance(asked, (int, float)) and not isinstance(asked, bool) and asked > cap:
        args["cap_ms"] = cap
    return {"check": "budget_respected", "args": args}


# ---------------------------------------------------------------------------------------------- relations
# Comparisons across two scenarios' answers, evaluated per arm by the comparator (lane S2b's compare.py); their
# semantics are stated in README.md, "Relations"
RELATION_TYPES = ("prefix_equals", "zero_equals_absent", "subset", "same_answer")


def relation(ctx: "Ctx", rtype: str, a: str, b: str, title: str, **args: Any) -> dict:
    """Record a relation between scenarios |a| and |b| (ids) on the context; the id is REL-<type>-nnnn."""
    if rtype not in RELATION_TYPES:
        raise ValueError("relation type %r" % rtype)
    n = sum(1 for r in ctx.relations if r["type"] == rtype) + 1
    rel = {"id": "REL-%s-%04d" % (rtype, n), "type": rtype, "a": a, "b": b, "title": title, "args": args}
    ctx.relations.append(rel)
    return rel


# ---------------------------------------------------------------------------------------------- scenarios
def scenario(sid: str, area: str, title: str, route: str, request: Any, truth: Optional[dict],
             expect: List[dict], timing_class: str = "fast", tags: Iterable[str] = (), repeat: int = 1,
             mcp: Optional[dict] = None, **extra: Any) -> dict:
    """One scenario of DESIGN §8.1. The MCP block is derived from the route and request unless given."""
    if mcp is None:
        if route == "POST /pattern":
            mcp = mcp_pattern(request)
        elif route == "POST /resolve":
            mcp = mcp_resolve(request)
        elif route == "POST /traverse":
            mcp = mcp_traverse(request)
        elif route == "GET /traverse/capabilities":
            mcp = {"tool": "traverse_capabilities", "arguments": {"databases": [DB]},
                   "by_design": ["capabilities.traverse_capabilities_document"]}
        elif route == "GET /capabilities":
            mcp = {"tool": "catalog", "arguments": {}, "by_design": ["capabilities.server_capabilities_document"]}
        elif route == "GET /pattern/capabilities":
            mcp = {"not_implemented": "the raw /pattern/capabilities document",
                   "by_design": ["capabilities.pattern_capabilities_document"]}
        else:
            mcp = {"not_implemented": "route %s" % route, "by_design": ["capabilities.not_expressible.protocol"]}
    if "body_raw" in extra or "headers" in extra or route.split(" ")[0] not in ("POST", "GET"):
        mcp = {"not_implemented": "protocol (raw body, headers or method)",
               "by_design": ["capabilities.not_expressible.protocol"]}
    # an exact scenario's source contexts are at level 0 (the MCP runner reads the level and substitutions)
    exact = not (isinstance(request, dict) and request.get("max_mismatches"))
    for e in expect:
        if e["check"] == "contains_source_context" and exact and "level" not in e["args"]:
            e["args"]["level"] = 0
    # every scenario that names a time budget and allows a 200 checks the elapsed time against it (a FINDING)
    if not any(e["check"] == "budget_respected" for e in expect):
        b = budget_expectation(route, request, expect)
        if b is not None:
            expect = list(expect) + [b]
    sc = {"id": sid, "area": area, "title": title, "route": route, "request": request, "mcp": mcp,
          "truth": truth or {}, "expect": expect, "timing_class": timing_class, "repeat": repeat,
          "tags": sorted(set(tags))}
    sc.update(extra)
    return sc


def source_contexts(t: dict, k: int = K) -> List[Tuple[str, int]]:
    """(kmer, offset) of every record k-mer holding the truth window, within the truth's flank."""
    out = []
    L = t["length"]
    rel = t["start"] - t["flank_start"]
    for p in range(0, k - L + 1):
        a = rel - p
        if a < 0 or a + k > len(t["flank"]):
            continue
        out.append((t["flank"][a:a + k], p))
    return out
