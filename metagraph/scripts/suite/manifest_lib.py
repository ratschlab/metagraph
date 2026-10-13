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
import random
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import seqtools as st
from truth import K, Record

DB = "refseq33m-experimental"

# Every pattern request names these (async finding 1: the MCP defaults differ: mode count, max_contexts 1,000)
PATTERN_TIME_BUDGET_MS = 60000
PATTERN_MAX_CONTEXTS = 10000
PATTERN_MAX_PREDICATE_CONTEXTS = 10000   # the MCP's default and the staging host's cap
# Explicit-label resolves name the deadline and the query bound the service sends (its max_time_ms 20,000 on
# staging and its own max_query_bp 2,000 / 500 for a discovery), so that both arms send the same body
RESOLVE_TIME_BUDGET_MS = 20000
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
# overrides traverse_fetch can express (DESIGN §4; async's traverse.not_expressible.strategy_keys)
TRAVERSE_OVERRIDABLE = {"bounds", "direction", "exhaustive", "labels.max_labels_per_node", "labels.seed_label_kind"}

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
        the window as the record holds it on its forward strand, and a flank to recompute the source contexts."""
        seq = rec.seq
        a = max(0, s - flank)
        b = min(len(seq), s + L + flank)
        return {
            "record": rec.accession, "genome": rec.key, "taxid": rec.column, "seq_id": rec.seq_id,
            "record_length": rec.length, "start": rec.offset + s, "length": L, "strand": strand,
            "window": seq[s:s + L], "flank": seq[a:b], "flank_start": rec.offset + a,
            "fetched_span": [rec.offset, rec.offset + len(seq)],
        }

    def flank_for(self, rec: Record, t: dict, extra: int) -> dict:
        """The same truth with a flank of |extra| bases each side (traversal walks)."""
        s = t["start"] - rec.offset
        out = dict(t)
        a = max(0, s - extra)
        b = min(len(rec.seq), s + t["length"] + extra)
        out["flank"] = rec.seq[a:b]
        out["flank_start"] = rec.offset + a
        return out


# ---------------------------------------------------------------------------------------------- requests
def pattern_request(patterns: List[dict], mode: str = "all_or_count", max_contexts: int = PATTERN_MAX_CONTEXTS,
                    time_budget_ms: int = PATTERN_TIME_BUDGET_MS, **fields: Any) -> dict:
    """A /pattern body that names mode, max_contexts and time_budget_ms (async finding 1)."""
    req = {"patterns": patterns, "mode": mode, "max_contexts": max_contexts, "time_budget_ms": time_budget_ms}
    for k, v in fields.items():
        if v is not None:
            req[k] = v
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
    """A /traverse body: the MCP preset's merged strategy, then the overrides, then |extra| (fields the MCP
    cannot express), with output.detail."""
    base = WALK_PRESET if preset == "walk" else TRIE_PRESET
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
    else:
        if labels == "predicate_only":
            ids.append("pattern.predicate_only_without_predicate")
        if "predicate_strands" in req or "predicate_scope" in req:
            ids.append("pattern.predicate_options_without_predicate")
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
        args["high_throughput"] = req["in_ram"]
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
        if len(labels) != len(set(labels)):
            ids.append("resolve.duplicate_labels")
        if len(set(labels)) > 1000:
            ids.append("resolve.too_many_labels")
        if any(isinstance(l, str) and l.strip() == "" for l in labels):
            ids.append("resolve.blank_label")
    if isinstance(disc, dict) and isinstance(disc.get("max_labels"), int) and not 1 <= disc["max_labels"] <= 100:
        ids.append("resolve.discover_limit_above_100")
    tb = (req.get("bounds") or {}).get("time_budget_ms")
    if isinstance(tb, (int, float)) and not isinstance(tb, bool):
        if tb > 20000:
            ids.append("resolve.time_budget_above_host_max")
        elif tb < 500:
            ids.append("resolve.time_budget_below_minimum")
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
        preset_name = "trace"
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
