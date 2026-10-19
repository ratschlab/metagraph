"""W5, the suite's independent oracle of the graphlet processing tools (manifest v8, DESIGN-v8 §2), phase A.

    answer(tool, args_resolved, full_answer, mgt_text, context=None)
        -> {"answer": {...} | None,
            "oracle": {"checked": {pointer: clause}, "consistency_only": [pattern, ...], "inputs": [...],
                       "declared"?: str, "note"?: str}}

The answer a library tool should give on the graphlet whose MGT text is |mgt_text| (W1's document: the server body,
or the standalone text with its J line) and whose `detail: full` twin is |full_answer| (the direct-local response,
or one result of it), for the arguments W1 sent (FORMAT-v8 §2: W5 takes W1's args_resolved). `answer` holds only
what W5 derives; `checked` names every value of it by its RFC 6901 pointer with the SPEC/DESIGN clause it rests on;
`consistency_only` names, as pointer patterns (`*` = any list index or object key), the fields of the tool's answer
W5 does not check: those are judged only by the agreement of W1..W4 (FORMAT-v8 §4). Two extensions (FORMAT
questions): `unordered` lists checked pointers whose list is compared as a multiset (a tool list no clause orders),
`absent` maps pointers that must not exist in the other answer to their clause. `differences(record, other)` is the
comparison the compare stage and the offline layer run; `export_view` turns an export's receipt and file into the
shape W5 answers.

`declared` says why a record checks nothing: `not_implemented_phase_b` (graphlet_support, graphlet_claims,
graphlet_labels name=: DESIGN-v8 §8 phase B), `properties_only` (graphlet_compare: the compare_* functions below),
`not_applicable` (backend and store tools), `view_consistency_only` (a tool on a subtrie view), `twin_mismatch`
(the full twin is not the graphlet's walk: twin_difference, the DESIGN-v8 §2 twin pre-check; nothing W5 says is
attributable then) and `expected_error` (the arguments name no answer, e.g. a walk id out of range: `differences`
then requires an error, of any code: the codes are the library's vocabulary, not the SPEC's).

`context` (all optional): seed_index (0); paging ("all", the default: lists are the whole list, every page
concatenated; otherwise rows are checked up to `rows`, the reference way's row count, and are consistency-only
without it); max_bytes (2048) and sequence_max_bytes (16384), the configured ceilings (DESIGN-v8 §3: from the
capabilities' limits block); view (the op runs on a subtrie handle); strategy (the graphlet envelope's normalized
strategy, for traverse_continue; default the full twin's).

Independence (DESIGN-v8 §2, REVIEW-DESIGN-v8 D-R5): this module imports the standard library and the suite's own
`checks` (its MGT reader, stdlib only) and nothing of `metagraph`; it cites SPEC-core
(docs/SPEC-labeled-traversal-core.md, "SPEC") and DESIGN (docs/DESIGN-traverse-graphlet.md) only. Every value comes,
by preference, from the engine's `full` twin (segments, paths with end_labels and continuation, runs, splits,
label_summary, limitations: computed by the C++ walker) and the MGT tokens (H, S); a field whose only definition is
a tool's own presentation is cited "INVENTORY §2 (presentation)" (list caps, tie orders, string forms), and a field
with no definition outside the library, or only DESIGN pseudocode the twin does not compute, is consistency-only.

Values a clock measured (a time limitation's observed, a time resource stop's amounts) differ between the twin's
execution and the graphlet's (SPEC §6.8) and are never checked.

Coverage (phase A; [S] SPEC, [D] DESIGN, [P] presentation per INVENTORY §2; c/o = consistency-only):
  every answer's evidence: complete_to_bp, scope [S §7.4], support, reconverge [S §7.5.2 H], outcome [S §7.0]
      | c/o: limitations (the tool's digest, finding F1), exact, coordinates, view.of
  graphlet_summary summary: seed.length_bp [S §7.1], seed.mode [S §7.4 label_mode], outcome [S §7.0],
      index.ns, index.verifiable [S §7.5.2 H, §7.5.7], arms.*: status, complete_to_bp, scope [S §7.4], walks,
      segments, splits, merges, bases, max_bp [S §7.5.5 counts, §7.5.2 A] | c/o: seed.seed_id, seed.labels,
      top_labels, caveats, cut flags
  graphlet_summary limitations: outcome, resource_stop [S §7.0], rows (each entry, arm) and total [S §7.0, §7.5.2 K
      order] | c/o: rows.*.informational; rows with arm= on a result with seed-level entries
  graphlet_walks: arm, total, rows.*.walk (in the rank's order: support [D §5], length, id), length_bp, end
      [S §7.5.6, §7.1; string P], n_alive [S §7.5.6], n_full, labels_full (+_truncated) [D §5.1 via
      end_labels.route_bp and runs.from_bp], sequence (spell tail/full) [S §7.5.4], label filter (alive at the leaf)
      and route_consistent (filtered.merge_entered: entered through a non-first parent) [D §6, §5.1; finding F2]
      | c/o: complete, beyond_certified_bp, coordinates, hint, filtered.not_from_seed; rank=loss order (no key in
      DESIGN); annotate: n_alive, n_full, labels_full, rank=support order, label filters
  graphlet_walk: arm, walk, length_bp, total, rows.*: segment, from_bp, to_bp, labels_in(_total), labels_out(_total)
      [S §7.1, §7.4], in, out [S §7.1; first 8 P], merge_of (absent on one parent) [S §7.5.6], ends (a set, absent
      without one) [S §7.5.6 label_end events; finding F5], continuation {length_bp, loss_used, labels} (absent
      without one) [S §7.1] | c/o: rows.*.exact; ends beyond 8; labels_out_total on cut lists (constrain)
  graphlet_labels (table): total, rows (name, ref, arm, direct_bp, reach_bp, reentries, runs) [S §7.1
      label_summary], order (-rank, ref, arm) [D §6; ties P], min_direct_bp, at_bp (constrain: a run holds the
      base) [S §7.1 runs] | c/o: at_bp in annotate
  graphlet_splits: total, rows (at_bp, segment, kind, labels_before, branches.char, branches.labels =
      labels_distinct) in (at_bp, first child) order [S §7.4, §7.5.6], min_labels_before
  graphlet_sequence: arm, orientation, with_seed, from, length, sequence [S §7.5.4], to (an uncut slice) | c/o:
      ceiling_bytes, truncated/next_from/to of a cut slice, left-arm walk orientation with the seed
  graphlet_export: fasta records, ids and sequences [D §5 to_fasta, S §7.5.4]; gfa P count and the walks the paths
      spell [D §9]; json = the normalised full result per key [S §7.5.6, D §9 T37]; json what=walk chains;
      mgt body digest [S §7.5.2, FORMAT §3] | c/o: header tags, gfa S/L lines, records of gfa/json/mgt, what=compare
  graphlet_subtrie: view selectors, mode, arm, summary labels, per-arm complete_to_bp, scope, qualified [D §5
      subgraph] | c/o: the view's walk and segment counts; every tool on the view handle
  graphlet_save / graphlet_load: the saved text's body digest and its graphlet_arms (round trip) [S §7.5.2];
      load's summary as graphlet_summary's
  traverse_continue (execute=False): request.seeds (sequence; labels, absent in annotate mode [S §6.9, §5;
      finding F3]), the reduced loss budget and branch allowance (the block where the allowance is > 0, finding F6),
      parent arm/walk/overlap_bp [S §7.1 continuation, D §5.3] | c/o: labels.extra, notes, the rest of the strategy
      (direction, detail); execute=True (backend)
  graphlet_compare: properties only (compare_self_properties, compare_swap_properties, compare_walks_properties,
      compare_prefix_subset_properties) [D §5, §5.2.5, §3.1, S §6.10, §7.5.7]
  phase B (declared): graphlet_support, graphlet_claims, graphlet_labels name=
"""

from __future__ import annotations

import copy
import hashlib
import json
import re
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import checks

ARMS = ("left", "right")
_ARM_NAMES = {"left": "left", "right": "right", "l": "left", "r": "right"}
_SUPPORT = {"k": "kmer", "t": "trace"}            # SPEC §7.5.2 H <support k|t>
_RECONVERGE = {"m": "merge", "k": "keep"}         # SPEC §7.5.2 H <reconverge m|k>
LIST_CAP = 8      # the label lists a row shows (labels_full, in, out, ends): INVENTORY §2.6/§2.7 (presentation)
_SEQ_MARGIN = 2048   # bytes of a sequence page that are not the slice (evidence, handle, keys): a slice shorter
                     # than the ceiling by this much is answered uncut

PHASE_B = ("graphlet_support", "graphlet_claims")
NOT_APPLICABLE = ("traverse_capabilities", "traverse_resolve", "traverse_fetch", "graphlet_list", "graphlet_free")

# clauses cited by more than one tool
C_PATHS = "SPEC §7.5.6 paths (leaf ordinal in segment order, chain through parents[0]): the twin's paths[].id"
C_LENGTH = "SPEC §7.5.6 paths: length_bp = the leaf's from_bp + length_bp (the twin's paths[].length_bp)"
C_SEG = "SPEC §7.1 segments[] (id, parents, from_bp, length_bp); to_bp = from_bp + length_bp (SPEC §7.5.4)"
C_LABEL = "DESIGN §5 Graphlet.label: {name, ref}, ref h:<column>:<seq_id> | c:<column> (SPEC §7.1 label_dict)"
C_LABEL_LIST = C_LABEL + "; ascending label id (SPEC §7.5.1 RANGES, §7.1 lists); first 8: INVENTORY §2 (presentation)"
C_CTB = "SPEC §7.4 complete_to_bp per arm (§6.10)"
C_SCOPE = "SPEC §7.4 completeness_scope per arm"
C_OUTCOME = "SPEC §7.0 outcome (§7.5.2 O)"
C_ARM_ECHO = "the arm asked (DESIGN §6), short forms l/r named left/right (INVENTORY §1.5)"
C_SPELL = "SPEC §7.5.4 natural orientation (right flank root -> leaf, left flank reversed) from the twin's sequences"


C_BAD_ARM = "library-defined: an arm not retrieved, or none of two, is bad_arm (INVENTORY §1.5)"
C_NO_BASES = "SPEC §7.1 sequences: false omits the strings; the code no_bases is library-defined (INVENTORY §2.12)"
C_BAD_ARG = "library-defined: an argument outside the tool's domain is bad_argument (INVENTORY §1.1)"


class _Refused(Exception):
    """The arguments name no answer (an arm not retrieved, a walk out of range, an unknown label, ...): the tool must
    answer an error (DESIGN §6: errors are structured results). |code| is the error code W5 expects (None: any code),
    |clause| where the refusal and its code come from; most codes are the library's vocabulary ("library-defined")."""

    def __init__(self, why: str, code: Optional[str] = "bad_argument", clause: str = C_BAD_ARG):
        super().__init__(why)
        self.code, self.clause = code, clause


def _int_arg(args: dict, key: str, default: Any, low: int) -> Any:
    """An integer argument (never a boolean) of at least |low|, or its default; else a refusal (library-defined
    argument checks, INVENTORY §1.1, §2.6, §2.9, §2.10)."""
    v = args.get(key, default)
    if v is default and default is None:
        return v
    if not isinstance(v, int) or isinstance(v, bool) or v < low:
        raise _Refused("%s %r: an integer >= %d" % (key, v, low))
    return v


def _bool_arg(args: dict, key: str, default: bool) -> bool:
    v = args.get(key, default)
    if not isinstance(v, bool):
        raise _Refused("%s %r: a boolean" % (key, v))
    return v


# ---------------------------------------------------------------------------------------------- the document
_O_CODES = ({"c": "complete", "p": "partial", "f": "failed"}, {"c": "complete", "x": "cut"},
            {"c": "complete", "l": "lower_bound", "q": "qualified"}, {"i": "inline", "s": "spooled", "p": "paged"})


def _mgt_result(mgt_text: str) -> dict:
    """A result read from the MGT text alone, as far as the compare properties need it (a body-only document has no
    twin): checks.graphlet_arms' arms with complete_to_bp and completeness_scope from the A records, outcome from O
    (SPEC §7.5.2)."""
    arms = checks.graphlet_arms(mgt_text)
    res: dict = {"arms": arms}
    for line in mgt_text.split("\n"):
        f = line.split(" ")
        if f[0] == "O" and len(f) >= 5:
            res["outcome"] = dict(zip(("walks", "branch_diagnostics", "label_evidence", "delivery"),
                                      (c.get(x, x) for c, x in zip(_O_CODES, f[1:5]))))
        elif f[0] == "A" and len(f) >= 5:
            a = arms.get({"l": "left", "r": "right"}.get(f[1], f[1]))
            if a is not None:
                a["complete_to_bp"] = int(f[3])
                a["completeness_scope"] = {"p": "per_path", "u": "united_history"}.get(f[4], f[4])
    return res


class _Doc:
    """The twin result, its envelope and the MGT tokens W5 reads (H, S), with the label dictionary. Without a twin
    (|full_answer| None, the compare properties on a body-only document) the result is read from the text."""

    def __init__(self, full_answer: Optional[dict], mgt_text: str, seed_index: int = 0):
        if full_answer is None:
            self.envelope = {}
            self.result = _mgt_result(mgt_text)
        elif isinstance(full_answer, dict) and isinstance(full_answer.get("results"), list):
            self.envelope = full_answer
            self.result = full_answer["results"][seed_index]
        else:
            self.envelope = {}
            self.result = full_answer
        self.mgt = mgt_text
        self.H: List[str] = []
        self.S: List[str] = []
        for line in mgt_text.split("\n"):
            if line.startswith("H ") and not self.H:
                self.H = line.split(" ")
            elif line.startswith("S "):
                self.S = line.split(" ", 5)
                break
        self.arms: Dict[str, dict] = {s: a for s, a in (self.result.get("arms") or {}).items() if s in ARMS}
        self.sides = [s for s in ARMS if s in self.arms]
        self.mode = self.result.get("label_mode") or {"c": "constrain", "a": "annotate"}.get(self._h(6), self._h(6))
        self.labels = self._label_list()

    def _h(self, i: int) -> Optional[str]:
        """H field |i| (SPEC §7.5.2: H mgt 1 k regime alphabet mode support reconverge cap continuation_bp seed_index
        walk index_ns index_fp index_meta_fp), None for `*` or a missing field."""
        v = self.H[i] if i < len(self.H) else None
        return None if v in (None, "*") else v

    @property
    def seed(self) -> str:
        """The seed as validated, request orientation (SPEC §7.5.2 S, its last field)."""
        return self.S[5] if len(self.S) > 5 else ""

    def _label_list(self) -> List[dict]:
        """[{name, ref}] per label id: the twin's label_dict (SPEC §7.1: column, and seq_id for header labels), else
        the L records (SPEC §7.5.2)."""
        ld = self.result.get("label_dict") or []
        if not ld or any("column" not in e for e in ld):
            ld = checks.graphlet_label_dict(self.mgt)
        out = []
        for e in ld:
            ref = "h:%s:%s" % (e.get("column"), e.get("seq_id")) if e.get("kind") == "header" else "c:%s" % e.get("column")
            out.append({"name": e["name"], "ref": ref})
        return out

    def label(self, i: int) -> dict:
        return dict(self.labels[i])

    def resolve(self, sel: Any) -> int:
        """A label selector -> its id (DESIGN §5 Graphlet.label: tagged {ref} | {name} | {id}; a bare string resolves
        when exactly one label matches it as a ref or as a name)."""
        c_lab = "DESIGN §5 Graphlet.label (AmbiguousLabel), §6 unknown_label; the codes are library-defined " \
                "(INVENTORY §1.1)"
        hits: List[int] = []
        if isinstance(sel, dict) and len(sel) == 1:
            (key, val), = sel.items()
            if key == "id" and isinstance(val, int) and not isinstance(val, bool) and 0 <= val < len(self.labels):
                return val
            if key in ("ref", "name"):
                hits = [i for i, x in enumerate(self.labels) if x[key] == val]
            elif key != "id":
                raise _Refused("label selector %r" % (sel,), None, c_lab)
        elif isinstance(sel, str):
            hits = sorted({i for i, x in enumerate(self.labels) if sel in (x["ref"], x["name"])})
        else:
            raise _Refused("label selector %r" % (sel,), None, c_lab)
        if len(hits) == 1:
            return hits[0]
        raise _Refused("label selector %r: %s" % (sel, "%d labels" % len(hits) if hits else "no label"),
                       "ambiguous_label" if hits else "unknown_label", c_lab)

    def arm(self, arg: Any, required: bool) -> Optional[str]:
        """The arm an argument names (left/right/l/r); None names the only arm of a one-arm retrieval, or every arm
        where the tool allows it (INVENTORY §1.5)."""
        if arg is None:
            if len(self.sides) == 1:
                return self.sides[0]
            if required:
                raise _Refused("both arms are present: the tool needs one", "bad_arm", C_BAD_ARM)
            return None
        side = _ARM_NAMES.get(arg) if isinstance(arg, str) else None
        if side is None or side not in self.arms:
            raise _Refused("arm %r is not an arm of this retrieval" % (arg,), "bad_arm", C_BAD_ARM)
        return side

    def flanks(self, side: str) -> Optional[List[str]]:
        """Per path id, the flank in natural orientation (SPEC §7.1, §7.5.4), or None without sequences."""
        try:
            return [x for _, x in checks._arm_flanks(self.arms[side], side)]
        except checks.Missing:
            return None

    def molecule(self, side: str, flank: str) -> str:
        """The walk with the seed, natural orientation: natural(left) + seed or seed + natural(right) (SPEC §7.5.4)."""
        return flank + self.seed if side == "left" else self.seed + flank


# ---------------------------------------------------------------------------------------------- the record
def _esc(key: Any) -> str:
    return str(key).replace("~", "~0").replace("/", "~1")


def _matches(pattern: str, pointer: str) -> bool:
    p, q = pattern.split("/"), pointer.split("/")
    return len(p) == len(q) and all(a == "*" or a == b for a, b in zip(p, q))


def _record(ans: Optional[dict], clauses: Dict[str, str], consistency: Sequence[str],
            inputs: Sequence[str] = ("full", "mgt"), unordered: Sequence[str] = (), absent: Sequence[str] = (),
            **extra: Any) -> dict:
    """The W5 record: every value of |ans| claimed by the first clause pattern that matches its pointer (a pattern
    naming an object or a list claims it whole). A value no clause claims is a bug of this module: W5 answers only
    what it checks (FORMAT-v8 §4). |unordered|: patterns of checked lists compared as multisets (a tool list whose
    order no clause fixes); |absent|: pointers (with their clause) that must not exist in the other answer, e.g. a
    seed's labels in annotate mode (FORMAT questions Q2, Q4)."""
    checked: Dict[str, str] = {}

    def visit(node: Any, ptr: str) -> None:
        for pat, clause in clauses.items():
            if _matches(pat, ptr):
                checked[ptr] = clause
                return
        if isinstance(node, dict) and node:
            for k, v in node.items():
                visit(v, ptr + "/" + _esc(k))
        elif isinstance(node, list) and node:
            for i, v in enumerate(node):
                visit(v, ptr + "/" + str(i))
        elif not (isinstance(node, (dict, list)) and not node and ptr):
            raise AssertionError("graphlet_oracle: no clause claims %s" % (ptr or "/"))

    if ans is not None:
        for k, v in ans.items():
            visit(v, "/" + _esc(k))
    oracle: dict = {"checked": checked, "consistency_only": sorted(set(consistency)), "inputs": list(inputs)}
    unord = sorted(p for p in checked if any(_matches(u, p) for u in unordered))
    if unord:
        oracle["unordered"] = unord
    if absent:
        oracle["absent"] = dict(absent)
    oracle.update({k: v for k, v in extra.items() if v is not None or k == "expected_error"})
    return {"answer": ans, "oracle": oracle}


def _declared(why: str, note: Optional[str] = None, inputs: Sequence[str] = ("full", "mgt")) -> dict:
    return _record(None, {}, ["/*"], inputs, declared=why, note=note)


def resolve_pointer(doc: Any, pointer: str) -> Any:
    """RFC 6901: the value at |pointer| in |doc|; KeyError when absent."""
    cur = doc
    if pointer in ("", "/"):
        return cur
    for raw in pointer.split("/")[1:]:
        tok = raw.replace("~1", "/").replace("~0", "~")
        if isinstance(cur, dict):
            if tok not in cur:
                raise KeyError(pointer)
            cur = cur[tok]
        elif isinstance(cur, list):
            if not tok.isdigit() or int(tok) >= len(cur):
                raise KeyError(pointer)
            cur = cur[int(tok)]
        else:
            raise KeyError(pointer)
    return cur


def same(a: Any, b: Any) -> bool:
    """JSON value equality as the suite compares parsed values: booleans only with booleans, numbers by value
    (an integral float equals its integer), objects by keys, lists in order."""
    if isinstance(a, bool) or isinstance(b, bool):
        return isinstance(a, bool) and isinstance(b, bool) and a == b
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        return float(a) == float(b)
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(same(a[k], b[k]) for k in a)
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(same(x, y) for x, y in zip(a, b))
    return type(a) is type(b) and a == b


def error_code(other: Any) -> Optional[str]:
    """The error code of a tool answer ({"error": "<code>", ...}) or of a per-call record's error ({"code", "fields"},
    FORMAT v2 B2); None when |other| is no error."""
    if not isinstance(other, dict):
        return None
    e = other.get("error")
    if isinstance(e, str):
        return e
    if isinstance(e, dict) and isinstance(e.get("code"), str):
        return e["code"]
    return None


def _canon(v: Any) -> str:
    return json.dumps(v, sort_keys=True, separators=(",", ":"))


def _cut_by_other(other: Any, ptr: str) -> bool:
    """Whether |other| names the field at |ptr| as cut to fit its page (DESIGN §6: a row larger than max_bytes comes
    alone with row_truncated and the fields cut named): on the object that holds the field, or, for a page of one
    row, on the page (FORMAT question Q5: an assembled answer keeps only the first page's marker)."""
    parent, key = ptr.rsplit("/", 1)
    try:
        holder = resolve_pointer(other, parent)
    except KeyError:
        return False
    if isinstance(holder, dict) and holder.get("row_truncated") and key in (holder.get("cut_fields") or ()):
        return True
    m = re.match(r"/rows/0/[^/]+\Z", ptr)
    return bool(m and isinstance(other, dict) and other.get("row_truncated") and len(other.get("rows") or []) == 1
                and key in (other.get("cut_fields") or ()))


def _is_cut_of(mine: Any, theirs: Any, unordered: bool) -> bool:
    """|theirs| a cut of |mine| in the suite's words (REVIEW D-R8): a string or a list a prefix, an unordered list a
    sub-multiset."""
    if isinstance(mine, str) and isinstance(theirs, str):
        return mine.startswith(theirs)
    if isinstance(mine, list) and isinstance(theirs, list):
        if unordered:
            pool = [_canon(x) for x in mine]
            for x in theirs:
                if _canon(x) not in pool:
                    return False
                pool.remove(_canon(x))
            return True
        return len(theirs) <= len(mine) and all(same(x, y) for x, y in zip(mine, theirs))
    return same(mine, theirs)


def differences(record: dict, other: Optional[dict]) -> List[dict]:
    """Where |other| (another way's assembled answer, service fields removed) departs from W5 on a checked
    pointer: [{pointer, clause, oracle, other | missing | present}]. Lists the record names `unordered` are compared
    as multisets; a field |other| names cut to fit its page must be a cut of W5's (_is_cut_of); a pointer the record
    names `absent` must not resolve in |other|. A record that declares `expected_error` requires an `error` key in
    |other|; any other declaration checks nothing."""
    o = record["oracle"]
    if o.get("declared") == "expected_error":
        code = error_code(other)
        if code is None:
            return [{"pointer": "/error", "clause": "DESIGN §6: errors are structured results (%s)" % o.get("note"),
                     "oracle": o.get("expected_error") or "an error", "missing": True}]
        if o.get("expected_error") and code != o["expected_error"]:
            return [{"pointer": "/error/code", "clause": o.get("note"), "oracle": o["expected_error"], "other": code}]
        return []
    out = []
    for ptr, clause in o["checked"].items():
        mine = resolve_pointer(record["answer"], ptr)
        try:
            theirs = resolve_pointer(other, ptr)
        except KeyError:
            out.append({"pointer": ptr, "clause": clause, "oracle": mine, "missing": True})
            continue
        unord = ptr in o.get("unordered", ())
        if _cut_by_other(other, ptr):
            if not _is_cut_of(mine, theirs, unord):
                out.append({"pointer": ptr, "clause": clause + "; cut to fit (DESIGN §6, REVIEW D-R8): a prefix",
                            "oracle": mine, "other": theirs})
            continue
        if unord and isinstance(mine, list) and isinstance(theirs, list):
            mine, theirs = sorted(mine, key=_canon), sorted(theirs, key=_canon)
        if not same(mine, theirs):
            out.append({"pointer": ptr, "clause": clause, "oracle": mine, "other": theirs})
    for ptr, clause in (o.get("absent") or {}).items():
        try:
            theirs = resolve_pointer(other, ptr)
        except KeyError:
            continue
        out.append({"pointer": ptr, "clause": clause, "oracle": "absent", "present": theirs})
    return out


# ---------------------------------------------------------------------------------------------- shared derivations
def twin_difference(full_answer: dict, mgt_text: str, seed_index: int = 0) -> Optional[str]:
    """The DESIGN-v8 §2 twin pre-check: checks.graphlet_arms of the graphlet body against the full twin's arms
    (status, segments: id, parents, from_bp, length_bp; paths: id, segments, length_bp, end_reasons, path_reason).
    None when they describe one walk; else the first difference."""
    res = full_answer["results"][seed_index] if isinstance(full_answer.get("results"), list) else full_answer
    try:
        arms = checks.graphlet_arms(mgt_text)
    except (checks.Missing, ValueError, KeyError, IndexError) as e:
        return "the MGT text does not read: %s" % e
    full_arms = {s: a for s, a in (res.get("arms") or {}).items() if s in ARMS}
    if set(full_arms) != set(arms):
        return "arms: full %s, graphlet %s" % (sorted(full_arms), sorted(arms))
    for side, fa in full_arms.items():
        ga = arms[side]
        if fa.get("status") != ga.get("status"):
            return "%s status: full %r, graphlet %r" % (side, fa.get("status"), ga.get("status"))
        fs = [(s["id"], s["parents"], s["from_bp"], s["length_bp"]) for s in fa.get("segments") or []]
        gs = [(s["id"], s["parents"], s["from_bp"], s["length_bp"]) for s in ga["segments"]]
        if fs != gs:
            i = next((i for i, (x, y) in enumerate(zip(fs, gs)) if x != y), min(len(fs), len(gs)))
            return "%s segment %d: full %s, graphlet %s (%d vs %d segments)" % (
                side, i, fs[i] if i < len(fs) else None, gs[i] if i < len(gs) else None, len(fs), len(gs))
        fp = [(p["id"], p["segments"], p["length_bp"], p.get("end_reasons") or {}, p.get("path_reason"))
              for p in fa.get("paths") or []]
        gp = [(p["id"], p["segments"], p["length_bp"], p["end_reasons"], p["path_reason"]) for p in ga["paths"]]
        if fp != gp:
            i = next((i for i, (x, y) in enumerate(zip(fp, gp)) if x != y), min(len(fp), len(gp)))
            return "%s path %d: full %s, graphlet %s (%d vs %d paths)" % (
                side, i, fp[i] if i < len(fp) else None, gp[i] if i < len(gp) else None, len(fp), len(gp))
    return None


def body_digest(mgt_text: str) -> str:
    """The body identity of FORMAT-v8 §3 (doc.body_digest): sha256 over the MGT records without the J line and the
    Z record (the J line is the library's envelope, the Z record counts it: SPEC §7.5.2), H's seed_index (the seed's
    position in the request, field 11) masked to `*`. A local equivalent of lane L's graphlet_doc.body_digest."""
    h = hashlib.sha256()
    for line in mgt_text.split("\n"):
        if not line or line.startswith("J ") or line.startswith("Z "):
            continue
        if line.startswith("H "):
            f = line.split(" ")
            if len(f) > 11:
                f[11] = "*"
            line = " ".join(f)
        h.update(line.encode("utf-8") + b"\n")
    return h.hexdigest()


def arms_digest(mgt_text: str) -> str:
    """sha256 of the canonical JSON of checks.graphlet_arms(mgt_text): the suite's reading of the whole graphlet
    (DAG, bases, entry sets, events, paths), for the save/load round trip."""
    text = json.dumps(checks.graphlet_arms(mgt_text), sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def normalise_full_result(result: dict) -> dict:
    """A `detail: full` result as two executions of one walk state it alike (the suite's normalisation, not the
    library's): `timing` removed everywhere and segment events sorted within equal at_bp, needed_budgets sorted
    (DESIGN §9 T37; SPEC §7.1 "events at one at_bp are unordered", §7.2 needed_budgets a histogram); growth
    label_ends null read as {} (DESIGN §2.5); annotation's direct_reads and keys_mapped dropped (SPEC §8.3: they
    count lookahead work, which differs between executions)."""
    r = checks._strip(copy.deepcopy(result), ())
    ann = r.get("annotation")
    if isinstance(ann, dict):
        for k in ("direct_reads", "keys_mapped"):
            ann.pop(k, None)
    for a in (r.get("arms") or {}).values():
        if not isinstance(a, dict):
            continue
        for s in a.get("segments") or []:
            if isinstance(s.get("events"), list):
                s["events"] = sorted(s["events"], key=lambda e: (e.get("at_bp", 0), json.dumps(e, sort_keys=True)))
        if isinstance(a.get("needed_budgets"), list):
            a["needed_budgets"] = sorted(a["needed_budgets"], key=lambda x: json.dumps(x, sort_keys=True))
        for b in a.get("growth") or []:
            if isinstance(b, dict) and b.get("label_ends") is None and "label_ends" in b:
                b["label_ends"] = {}
    return r


_TIME_KNOBS = ("bounds.time_budget_ms",)
_TIME_RESOURCES = ("time", "attempt_deadline", "cancelled")


def _time_measured(result: dict) -> List[str]:
    """Pointers (relative to a full result) of the values a clock measured: a time limitation's observed, a time
    resource stop's amounts and message. Two executions of one walk measure different times (SPEC §6.8: a stop by
    time is not deterministic; §7.2 timing), so W5 checks none of them against its twin."""
    out = []
    for j, e in enumerate(result.get("limitations") or []):
        if e.get("knob") in _TIME_KNOBS:
            out.append("/limitations/%d/observed" % j)
    for s, a in (result.get("arms") or {}).items():
        for j, e in enumerate(a.get("limitations") or [] if isinstance(a, dict) else []):
            if e.get("knob") in _TIME_KNOBS:
                out.append("/arms/%s/limitations/%d/observed" % (_esc(s), j))
    rs = result.get("resource_stop")
    if isinstance(rs, dict) and rs.get("resource") in _TIME_RESOURCES:
        out += ["/resource_stop/%s" % k for k in ("used", "remaining", "effective", "message") if k in rs]
    return out


def _evidence(d: _Doc, sides: Sequence[str]) -> Tuple[dict, Dict[str, str], List[str]]:
    """The evidence block of every local answer (DESIGN §6: complete_to_bp, exact, support, reconverge, scope,
    outcome, limitations), for the arms the answer covers."""
    ev: dict = {"complete_to_bp": {s: d.arms[s]["complete_to_bp"] for s in sides},
                "scope": {s: d.arms[s]["completeness_scope"] for s in sides},
                "outcome": d.result["outcome"]}
    if d._h(7) in _SUPPORT:
        ev["support"] = _SUPPORT[d._h(7)]
    if d._h(8) in _RECONVERGE:
        ev["reconverge"] = _RECONVERGE[d._h(8)]
    clauses = {"/evidence/complete_to_bp": C_CTB, "/evidence/scope": C_SCOPE, "/evidence/outcome": C_OUTCOME,
               "/evidence/support": "SPEC §7.5.2 H support k|t (kmer, trace)",
               "/evidence/reconverge": "SPEC §7.5.2 H reconverge m|k (merge, keep)"}
    # evidence.limitations is the tool's digest of the stated limitations (finding F1: it leaves out a scope entry
    # that observed 0 merges and adds the resource stop): DESIGN §6 names the key without defining it
    return ev, clauses, ["/evidence/exact", "/evidence/coordinates", "/evidence/view", "/evidence/qualified",
                         "/evidence/limitations"]


def _cut_rows(ans: dict, ctx: dict, consistency: List[str], key: str = "rows") -> None:
    """Whole lists unless the op read fewer pages: then the reference way's row count, or no rows at all."""
    if ctx.get("paging", "all") == "all":
        return
    if isinstance(ctx.get("rows"), int):
        ans[key] = ans[key][:ctx["rows"]]
    else:
        ans[key] = []
        consistency.append("/%s/*" % key)


def _arm_counts(a: dict) -> dict:
    """An arm's counts (SPEC §7.5.5 counts; §7.5.2 A's n_segments, n_leaves, n_splits, n_merges, bases)."""
    segs = a.get("segments") or []
    return {"walks": len(a.get("paths") or []), "segments": len(segs), "splits": len(a.get("splits") or []),
            "merges": sum(1 for s in segs if len(s.get("parents") or []) > 1),
            "bases": sum(s.get("length_bp", 0) for s in segs),
            "max_bp": max((p["length_bp"] for p in a.get("paths") or []), default=0)}


# ---------------------------------------------------------------------------------------------- the tools
def _summary_part(d: _Doc, sides: Sequence[str]) -> Tuple[dict, Dict[str, str], List[str]]:
    s = {"seed": {"length_bp": d.result["seed"]["length_bp"], "mode": d.mode}, "outcome": d.result["outcome"],
         "index": {"ns": d._h(13), "verifiable": d._h(14) is not None}, "arms": {}}
    for side in sides:
        a = d.arms[side]
        s["arms"][side] = dict({"status": a["status"], "complete_to_bp": a["complete_to_bp"],
                                "scope": a["completeness_scope"]}, **_arm_counts(a))
    clauses = {"/summary/seed/length_bp": "SPEC §7.1 seed.length_bp (§7.5.2 S)",
               "/summary/seed/mode": "SPEC §7.4 label_mode",
               "/summary/outcome": C_OUTCOME,
               "/summary/index/ns": "SPEC §7.5.2 H index_ns (null for *)",
               "/summary/index/verifiable": "SPEC §7.5.7: a join is verifiable only with index_fp (H)",
               "/summary/arms/*/status": "SPEC §7.4 status", "/summary/arms/*/complete_to_bp": C_CTB,
               "/summary/arms/*/scope": C_SCOPE,
               "/summary/arms/*/walks": "SPEC §7.5.5 counts.leaves (§7.5.2 A n_leaves): one walk per leaf",
               "/summary/arms/*/segments": "SPEC §7.5.5 counts.segments (A n_segments)",
               "/summary/arms/*/splits": "SPEC §7.5.5 counts.splits (A n_splits)",
               "/summary/arms/*/merges": "SPEC §7.5.5 counts.merges (A n_merges: segments with > 1 parent)",
               "/summary/arms/*/bases": "SPEC §7.5.5 counts.bases (A bases: the sum of G.length_bp)",
               "/summary/arms/*/max_bp": "SPEC §7.5.5 counts.max_bp (the longest path)"}
    cons = ["/summary/seed/seed_id", "/summary/seed/labels", "/summary/seed/derivation", "/summary/arms/*/top_labels",
            "/summary/caveats", "/summary/caveat_effects_cut", "/summary/caveats_truncated", "/summary/coordinates",
            "/parent", "/summary_stop", "/fields_cut"]
    return s, clauses, cons


def _graphlet_summary(d: _Doc, args: dict, ctx: dict) -> dict:
    side = d.arm(args.get("arm"), required=False)
    sides = [side] if side else d.sides
    detail = args.get("detail", "summary")
    if detail == "limitations":
        return _limitations(d, side, ctx)
    if detail != "summary":
        raise _Refused("detail %r: summary or limitations (DESIGN §6)" % (detail,))
    summ, clauses, cons = _summary_part(d, sides)
    ev, ecl, econ = _evidence(d, sides)
    clauses.update(ecl)
    return _record({"summary": summ, "evidence": ev}, clauses, cons + econ)


def _limitations(d: _Doc, side: Optional[str], ctx: dict) -> dict:
    """detail='limitations': every stated limitation in full, the seed's first, then the left arm's and the right
    arm's, each in its JSON order (SPEC §7.5.2 K), and the resource stop (SPEC §7.0)."""
    seed_level = list(d.result.get("limitations") or [])
    # (the entry's pointer in the result, the row): seed-level entries carry no arm
    src = [("/limitations/%d" % j, dict(e)) for j, e in enumerate(seed_level)] if side is None else []
    for s in d.sides:
        if side is None or s == side:
            src += [("/arms/%s/limitations/%d" % (s, j), dict(e, arm=s))
                    for j, e in enumerate(d.arms[s].get("limitations") or [])]
    rows = [r for _, r in src]
    timed = set(_time_measured(d.result))
    cons_t = []
    for i, (base, row) in enumerate(src):
        if base + "/observed" in timed:
            row.pop("observed", None)
            cons_t.append("/rows/%d/observed" % i)
    rs = d.result.get("resource_stop")
    c_rs = "SPEC §7.0 resource_stop (§7.5.2 Q), null without one"
    clauses = {"/outcome": C_OUTCOME, "/resource_stop": c_rs}
    if isinstance(rs, dict) and any(t.startswith("/resource_stop/") for t in timed):
        rs = {k: v for k, v in rs.items() if "/resource_stop/" + k not in timed}
        cons_t += [t for t in timed if t.startswith("/resource_stop/")]
        clauses = {"/outcome": C_OUTCOME, "/resource_stop/*": c_rs + "; a clock's amounts left out (SPEC §6.8)"}
    ans = {"outcome": d.result["outcome"], "resource_stop": rs, "total": len(rows), "rows": rows}
    clauses.update({"/total": "SPEC §7.0: every stated limitation (seed level and the arms asked)",
                    "/rows/*/arm": "SPEC §7.5.2 K arm (an arm's entry names its arm)",
                    "/rows/*/*": "SPEC §7.0 limitations: the stated entry (§7.5.2 K: the JSON entry itself)"})
    cons = ["/rows/*/informational", "/next_cursor"] + cons_t
    if side is not None and seed_level:
        # whether a seed-level entry belongs to an arm's listing is the tool's choice
        ans.pop("total")
        ans["rows"] = []
        cons += ["/total", "/rows/*"]
    _cut_rows(ans, ctx, cons)
    return _record(ans, clauses, cons)


def _evidence_from_leaf(e: dict, runs: List[dict]) -> int:
    """DESIGN §5.1 displayed-path evidence of a leaf label: max(the merge part, the run's from_bp), where the merge
    part of a leaf label equals its T.route_bp (the twin's end_labels[].route_bp; asserted by T37)."""
    run = runs[e["run"]] if isinstance(e.get("run"), int) and 0 <= e["run"] < len(runs) else {}
    return max(e.get("route_bp") or 0, run.get("from_bp") or 0)


def _end_text(p: dict) -> str:
    """A walk's end: its path_reason, else its labels' end reasons (sorted, comma-joined), else "semantic" (SPEC §7.1
    path_reason, end_reasons; §7.5.5 leaves_by_reason; the string form: INVENTORY §2.6)."""
    if p.get("path_reason"):
        return p["path_reason"]
    if p.get("end_reasons"):
        return ",".join(sorted(p["end_reasons"]))
    return "semantic"


def _graphlet_walks(d: _Doc, args: dict, ctx: dict) -> dict:
    side = d.arm(args.get("arm"), required=True)
    a = d.arms[side]
    rank = args.get("rank", "support")
    spell = args.get("spell", "none")
    if rank not in ("support", "length", "loss", "id") or spell not in ("none", "tail", "full"):
        raise _Refused("rank %r / spell %r" % (rank, spell))
    tail_bp = _int_arg(args, "tail_bp", 60, 0)
    min_bp = _int_arg(args, "min_bp", 0, 0)
    _int_arg(args, "n", 10, 1)
    rc = _bool_arg(args, "route_consistent", True)
    if _bool_arg(args, "coordinates", False) and not isinstance(d.result.get("coordinates"), dict):
        # SPEC §7.1: a retrieval without record coordinates states coordinates: null (with its reason) or none
        raise _Refused("coordinates=true on a retrieval without record coordinates (%s)" % (
            d.result.get("coordinates_reason") or "not requested"), "no_coordinates",
            "SPEC §7.1 record coordinates (null, coordinates_reason); the code no_coordinates is library-defined "
            "(INVENTORY §2.6)")
    lab = d.resolve(args["label"]) if args.get("label") is not None else None
    constrain = d.mode == "constrain"
    runs = a.get("runs") or []
    flanks = d.flanks(side) if spell != "none" else None
    cons = ["/rows/*/complete", "/rows/*/beyond_certified_bp", "/rows/*/coordinates", "/hint", "/next_cursor",
            "/filtered/not_from_seed"]
    rows, keys = [], []
    filtered = {"merge_entered": 0}
    for p in a["paths"]:
        if p["length_bp"] < min_bp:
            continue
        el = p.get("end_labels") or []
        full_ids = sorted(e["label"] for e in el if _evidence_from_leaf(e, runs) == 0)
        if lab is not None:
            if not constrain:
                break
            entry = next((e for e in el if e["label"] == lab), None)
            if entry is None:
                continue
            if rc and (entry.get("route_bp") or 0) > 0:
                # the label's kept entry came through a non-first parent at a merge on the displayed walk: the
                # displayed bases before it are not this lineage's (DESIGN §5.1; T.route_bp, finding F2)
                filtered["merge_entered"] += 1
                continue
        row: dict = {"walk": p["id"], "length_bp": p["length_bp"], "end": _end_text(p)}
        if constrain:
            row["n_alive"] = p.get("n_labels", len(el))
            row["n_full"] = len(full_ids)
            row["labels_full"] = [d.label(i) for i in full_ids[:LIST_CAP]]
            if len(full_ids) > LIST_CAP:
                row["labels_full_truncated"] = True
        if flanks is not None:
            f = flanks[p["id"]]
            row["sequence"] = f if spell == "full" else (f[:tail_bp] if side == "left" else f[max(0, len(f) - tail_bp):])
        rows.append(row)
        keys.append({"support": (-len(full_ids), -row.get("n_alive", 0), -p["length_bp"], p["id"]),
                     "length": (-p["length_bp"], p["id"]), "id": (p["id"],), "loss": (p["id"],)}[rank])
    rows = [r for _, r in sorted(zip(keys, rows), key=lambda kr: kr[0])]
    order = {"support": "DESIGN §5 walks(by='support'): (-|labels_full|, -n_alive, -length_bp, path_id)",
             "length": "DESIGN §5 walks(by='length'): longest first, ties by path_id",
             "id": "path id order", "loss": ""}[rank]
    ans: dict = {"arm": side, "total": len(rows), "rows": rows}
    clauses = {"/arm": C_ARM_ECHO,
               "/total": "DESIGN §5 walks(min_bp, labels): the arm's walks the filters keep (SPEC §7.5.6 paths; "
                         "label: alive at the leaf, end_labels)",
               "/rows/*/walk": C_PATHS + "; in the rank's order, " + order,
               "/rows/*/length_bp": C_LENGTH,
               "/rows/*/end": "SPEC §7.1 path_reason / end_reasons, §7.5.5 'semantic' (string form: INVENTORY §2.6)",
               "/rows/*/n_alive": "SPEC §7.5.6 end_labels: the labels alive at the leaf (the twin's paths[].n_labels)",
               "/rows/*/n_full": "DESIGN §5.1: leaf labels whose displayed evidence max(T.route_bp, run from_bp) is 0 "
                                 "(the twin's end_labels[].route_bp, runs[].from_bp)",
               "/rows/*/labels_full": "DESIGN §5.1 (as n_full); " + C_LABEL_LIST,
               "/rows/*/labels_full_truncated": "more than 8 full labels: INVENTORY §2.6 (presentation)",
               "/rows/*/sequence": C_SPELL + "; spell full = the flank, tail = its last tail_bp outward bases (DESIGN "
                                   "§6 spell=tail, tail_bp)",
               "/filtered/merge_entered": "DESIGN §6 'a filter never hides rows silently' (merge_entered), §5.1: "
                                          "with route_consistent, a walk the label reaches but entered through a "
                                          "non-first parent (the twin's end_labels[].route_bp > 0) is left out"}
    if lab is not None and rc and constrain and any(filtered.values()):
        ans["filtered"] = {k: v for k, v in filtered.items() if v}
    if rank == "loss":
        # DESIGN §5 names by='loss' but states no key: the order is the tool's
        cons.append("/rows/*")
        ans["rows"] = []
    if not constrain:
        # annotate: no lineages (SPEC §7.4: end_labels empty), displayed evidence is the witness rule of DESIGN
        # §5.2.1, which the twin does not compute
        cons += ["/rows/*/n_alive", "/rows/*/n_full", "/rows/*/labels_full", "/rows/*/labels_full_truncated"]
        if rank == "support":
            cons.append("/rows/*")
            ans["rows"] = []
        if lab is not None:
            ans = {"arm": side}
            cons += ["/total", "/rows/*", "/filtered"]
    if flanks is None and spell != "none":
        cons.append("/rows/*/sequence")
    if "rows" in ans:
        _cut_rows(ans, ctx, cons)
    ev, ecl, econ = _evidence(d, [side])
    ans["evidence"] = ev
    clauses.update(ecl)
    return _record(ans, clauses, cons + econ)


def _walk_rows(d: _Doc, side: str, p: dict, cap: Optional[int]) -> List[dict]:
    """One row per segment of the walk's first-parent chain (SPEC §7.5.6), the label lists cut to |cap|."""
    a = d.arms[side]
    segs = {s["id"]: s for s in a["segments"]}
    truncated_arm = (a.get("labels_per_node") or {}).get("nodes_truncated", 0) > 0
    rows = []
    for sid in p["segments"]:
        s = segs[sid]
        entry = s.get("labels") or []
        end = s.get("labels_at_end") or []
        row = {"segment": sid, "from_bp": s["from_bp"], "to_bp": s["from_bp"] + s["length_bp"],
               "labels_in": len(entry), "labels_in_total": s.get("labels_total", len(entry)),
               "labels_out": len(end), "in": [d.label(i) for i in entry[:cap]],
               "out": [d.label(i) for i in end[:cap]]}
        if d.mode == "annotate":
            ls = s.get("label_sets") or []
            row["labels_out_total"] = ls[-1]["labels_total"] if ls else row["labels_in_total"]
        elif not truncated_arm:
            row["labels_out_total"] = len(end)
        if len(s.get("parents") or []) > 1:
            row["merge_of"] = list(s["parents"])
        ends = [e for e in s.get("events") or [] if e.get("type") == "label_end"]
        if ends:
            # which `cap` of more ends a row shows depends on an order no clause fixes: such a row's ends are
            # consistency-only (_walk_cut_ends)
            row["ends"] = [dict(d.label(e["label"]), at_bp=e["at_bp"], why=e["reason"])
                           for e in sorted(ends, key=lambda e: (e["at_bp"], e["label"]))]
        rows.append(row)
    return rows


def _walk_row_extras(rows: List[dict], prefix: str, cap: Optional[int]) -> Tuple[List[str], Dict[str, str]]:
    """Per walk row: the consistency-only pointer of an `ends` list longer than |cap| (cut there), and the keys that
    must be absent: merge_of on a segment with one parent, ends on a segment with no label end."""
    cons, absent = [], {}
    for i, r in enumerate(rows):
        if cap is not None and len(r.get("ends") or []) > cap:
            del r["ends"]
            cons.append("%s/%d/ends" % (prefix, i))
        elif "ends" not in r:
            absent["%s/%d/ends" % (prefix, i)] = "SPEC §7.5.6: no R record ends a run on this segment (ends?)"
        if "merge_of" not in r:
            absent["%s/%d/merge_of" % (prefix, i)] = "SPEC §7.5.6: a segment with one parent merges nothing (merge_of?)"
    return cons, absent


_WALK_ROW_CLAUSES = {
    "segment": "SPEC §7.5.6 the walk's chain through parents[0] (the twin's paths[].segments)",
    "from_bp": C_SEG, "to_bp": C_SEG,
    "labels_in": "SPEC §7.1 segment.labels: the labels at the segment's first step (the recorded list)",
    "labels_in_total": "SPEC §7.4 segments[].labels_total (the true count)",
    "labels_out": "SPEC §7.1 labels_at_end: the labels at the segment's last step (§7.5.6 = G.end)",
    "labels_out_total": "SPEC §7.1 labels_at_end (constrain, no cut list) / §7.4 the last label_sets' labels_total "
                        "(annotate)",
    "in": "SPEC §7.1 segment.labels; " + C_LABEL_LIST, "out": "SPEC §7.1 labels_at_end; " + C_LABEL_LIST,
    "merge_of": "SPEC §7.5.6 reconverge: a segment with > 1 parent merges its parents (G parents)",
    "ends": "SPEC §7.5.6 label_end events (R with a code end, at to_bp, reason = qualifier text or enum), as a set "
            "(no clause orders this list); " + C_LABEL + "; checked when at most 8 (INVENTORY §2.7)",
}


def _graphlet_walk(d: _Doc, args: dict, ctx: dict) -> dict:
    side = d.arm(args.get("arm"), required=True)
    w = args.get("walk")
    paths = d.arms[side]["paths"]
    if not isinstance(w, int) or isinstance(w, bool) or not 0 <= w < len(paths):
        raise _Refused("walk %r of %d" % (w, len(paths)))
    p = paths[w]
    rows = _walk_rows(d, side, p, LIST_CAP)
    cons_r, absent = _walk_row_extras(rows, "/rows", LIST_CAP)
    ans: dict = {"arm": side, "walk": w, "length_bp": p["length_bp"], "total": len(rows), "rows": rows}
    clauses = {"/arm": C_ARM_ECHO, "/walk": C_PATHS, "/length_bp": C_LENGTH,
               "/total": "SPEC §7.5.6 the segments of the walk's chain"}
    clauses.update({"/rows/*/" + k: v for k, v in _WALK_ROW_CLAUSES.items()})
    c = p.get("continuation")
    if c:
        ans["continuation"] = {"length_bp": len(c.get("sequence") or ""), "loss_used": c["loss_used"],
                               "labels": len(c.get("labels") or [])}
        clauses["/continuation/length_bp"] = "SPEC §7.1 continuation.sequence: its length"
        clauses["/continuation/loss_used"] = "SPEC §7.1 continuation.loss_used (the smallest terminal loss)"
        clauses["/continuation/labels"] = "SPEC §7.1 continuation.labels: their number"
    else:
        absent["/continuation"] = "SPEC §7.1: only a leaf ended by the radius, a resource reason or beam_pruned " \
                                  "carries a continuation"
    cons = ["/rows/*/exact", "/next_cursor"] + cons_r
    if d.mode == "constrain" and (d.arms[side].get("labels_per_node") or {}).get("nodes_truncated", 0) > 0:
        cons.append("/rows/*/labels_out_total")
    _cut_rows(ans, ctx, cons)
    if len(ans["rows"]) < len(rows):
        absent = {k: v for k, v in absent.items() if not k.startswith("/rows/")}
    ev, ecl, econ = _evidence(d, [side])
    ans["evidence"] = ev
    clauses.update(ecl)
    return _record(ans, clauses, cons + econ, unordered=["/rows/*/ends"], absent=absent)


def _alive_at(a: dict, label: int, at_bp: int) -> bool:
    """A label is alive at outward base |at_bp| when one of its runs on the arm holds that base: from_bp <= at_bp <
    to_bp (SPEC §7.1 runs[], half-open intervals of outward bases)."""
    return any(r["label"] == label and r["from_bp"] <= at_bp < r["to_bp"] for r in a.get("runs") or [])


def _graphlet_labels(d: _Doc, args: dict, ctx: dict) -> dict:
    side = d.arm(args.get("arm"), required=False)
    sides = [side] if side else d.sides
    rank = args.get("rank", "direct_bp")
    if rank not in ("direct_bp", "reach_bp"):
        raise _Refused("rank %r" % (rank,))
    _int_arg(args, "n", 20, 1)
    min_direct = _int_arg(args, "min_direct_bp", 0, 0)
    at_bp = _int_arg(args, "at_bp", None, 0)
    if _bool_arg(args, "coordinates", False):
        raise _Refused("coordinates without name (INVENTORY §2.9: coordinates only with name=)")
    cons = ["/next_cursor"]
    rows = []
    for e in d.result.get("label_summary") or []:
        for s in sides:
            v = e.get(s)
            if not isinstance(v, dict) or v["direct_bp"] < min_direct:
                continue
            if at_bp is not None and d.mode == "constrain" and not _alive_at(d.arms[s], e["label"], at_bp):
                continue
            rows.append(dict(d.label(e["label"]), arm=s, direct_bp=v["direct_bp"], reach_bp=v["reach_bp"],
                             reentries=v["reentries"], runs=len(v.get("runs") or [])))
    rows.sort(key=lambda r: (-r[rank], r["ref"], r["arm"]))
    ans: dict = {"total": len(rows), "rows": rows}
    clauses = {"/total": "SPEC §7.1 label_summary: one row per label and arm the filters keep",
               "/rows/*/name": C_LABEL, "/rows/*/ref": C_LABEL,
               "/rows/*/arm": "SPEC §7.1 label_summary per arm; order (-rank, ref, arm): DESIGN §6 ranked by "
                              "direct_bp|reach_bp, ties by ref then arm: INVENTORY §2.9 (presentation)",
               "/rows/*/direct_bp": "SPEC §7.1 label_summary direct_bp (engine summarize)",
               "/rows/*/reach_bp": "SPEC §7.1 label_summary reach_bp (engine summarize)",
               "/rows/*/reentries": "SPEC §7.1 label_summary reentries",
               "/rows/*/runs": "SPEC §7.1 label_summary runs: their number"}
    if at_bp is not None and d.mode != "constrain":
        # annotate tracks no runs (SPEC §7.4): which labels are alive at a base is the tool's reading of P runs
        ans = {}
        cons += ["/total", "/rows/*"]
    else:
        if at_bp is not None:
            clauses["/total"] += "; at_bp: a run of the label holds the base (SPEC §7.1 runs[] [from_bp, to_bp))"
        _cut_rows(ans, ctx, cons)
    ev, ecl, econ = _evidence(d, sides)
    ans["evidence"] = ev
    clauses.update(ecl)
    return _record(ans, clauses, cons + econ)


def _graphlet_splits(d: _Doc, args: dict, ctx: dict) -> dict:
    side = d.arm(args.get("arm"), required=True)
    m = _int_arg(args, "min_labels_before", 2, 0)
    _int_arg(args, "n", 10, 1)
    rows = [{"at_bp": sp["at_bp"], "segment": sp["segment"], "kind": sp["kind"],
             "labels_before": sp["labels_before"],
             "branches": [{"char": b["char"], "labels": b["labels_distinct"]} for b in sp.get("branches") or []]}
            for sp in d.arms[side].get("splits") or [] if sp["labels_before"] >= m]
    ans: dict = {"arm": side, "total": len(rows), "rows": rows}
    clauses = {"/arm": C_ARM_ECHO,
               "/total": "SPEC §7.4 splits[] with labels_before >= min_labels_before (DESIGN §6)",
               "/rows/*/at_bp": "SPEC §7.4 splits[].at_bp, in the order of SPEC §7.5.6 ((at_bp, first child))",
               "/rows/*/segment": "SPEC §7.4 splits[].segment", "/rows/*/kind": "SPEC §7.5.6 kind from G.split",
               "/rows/*/labels_before": "SPEC §7.4 labels_before (§7.5.6)",
               "/rows/*/branches/*/char": "SPEC §7.4 branches[].char (the child's first base)",
               "/rows/*/branches/*/labels": "SPEC §7.4 branches[].labels_distinct (the true count)"}
    cons = ["/next_cursor"]
    _cut_rows(ans, ctx, cons)
    ev, ecl, econ = _evidence(d, [side])
    ans["evidence"] = ev
    clauses.update(ecl)
    return _record(ans, clauses, cons + econ)


def _graphlet_sequence(d: _Doc, args: dict, ctx: dict) -> dict:
    side = d.arm(args.get("arm"), required=True)
    walk, segment = args.get("walk"), args.get("segment")
    orientation = args.get("orientation", "natural")
    with_seed = bool(args.get("with_seed", False))
    start, end = args.get("start", 0), args.get("end")
    a = d.arms[side]
    if (walk is None) == (segment is None) or orientation not in ("natural", "walk") or \
            (with_seed and segment is not None) or not isinstance(start, int) or start < 0 or \
            (end is not None and (not isinstance(end, int) or end < start)):
        raise _Refused("walk %r / segment %r / orientation %r / with_seed %r / [%r, %r)" % (
            walk, segment, orientation, with_seed, start, end))
    if walk is not None:
        if not isinstance(walk, int) or not 0 <= walk < len(a["paths"]):
            raise _Refused("walk %r" % (walk,))
        flanks = d.flanks(side)
        if flanks is None:
            raise _Refused("no bases (sequences: false)", "no_bases", C_NO_BASES)
        natural = d.molecule(side, flanks[walk]) if with_seed else flanks[walk]
    else:
        if not isinstance(segment, int) or not 0 <= segment < len(a["segments"]):
            raise _Refused("segment %r" % (segment,))
        if "sequence" not in a["segments"][segment]:
            raise _Refused("no bases (sequences: false)", "no_bases", C_NO_BASES)
        natural = a["segments"][segment]["sequence"]
    # walking order is the stored order of G bases: the natural order on the right arm, reversed on the left (SPEC
    # §7.5.4); SPEC defines it for flanks, not for a molecule with the seed
    text = natural[::-1] if orientation == "walk" and side == "left" else natural
    if start > len(text):
        raise _Refused("start %d beyond the %d bases" % (start, len(text)))
    stop = len(text) if end is None else min(end, len(text))
    ans: dict = {"arm": side, "orientation": orientation, "with_seed": with_seed, "from": start, "length": len(text)}
    clauses = {"/arm": C_ARM_ECHO, "/orientation": "the orientation asked (DESIGN §6)",
               "/with_seed": "the with_seed asked (DESIGN §6)", "/from": "the slice's start (DESIGN §6 graphlet_sequence)",
               "/length": "the length of the spelled walk or segment (SPEC §7.5.4; " + C_SPELL + ")",
               "/sequence": C_SPELL + "; walk orientation = the G bases' walking order (SPEC §7.5.4); with_seed = "
                            "natural(left) + seed | seed + natural(right); the slice [start, end) (DESIGN §6)",
               "/to": "the end of an uncut slice: min(end, length)"}
    cons = ["/ceiling_bytes"]
    fits = stop - start <= ctx.get("sequence_max_bytes", 16384) - _SEQ_MARGIN
    if fits or ctx.get("paging", "all") == "all":
        ans["sequence"] = text[start:stop]
    else:
        cons.append("/sequence")
    if fits:
        ans["to"] = stop
    else:
        cons += ["/to", "/truncated", "/next_from"]
    if orientation == "walk" and side == "left" and with_seed:
        ans.pop("sequence", None)
        cons.append("/sequence")
    ev, ecl, econ = _evidence(d, [side])
    ans["evidence"] = ev
    clauses.update(ecl)
    return _record(ans, clauses, cons + econ)


# ---------------------------------------------------------------------------------------------- exports
def parse_fasta(text: str) -> List[dict]:
    """[{id, header, sequence}] of a FASTA text (lines joined, a record's id = the header's first token)."""
    out: List[dict] = []
    for line in text.split("\n"):
        if line.startswith(">"):
            out.append({"id": line[1:].split(" ", 1)[0], "header": line[1:], "sequence": []})
        elif out and line.strip():
            out[-1]["sequence"].append(line.strip())
    for r in out:
        r["sequence"] = "".join(r["sequence"])
    return out


def _cigar_overlap(c: str) -> Optional[int]:
    return int(c[:-1]) if c and c[-1] == "M" and c[:-1].isdigit() else (0 if c in ("0M",) else None)


def gfa_paths(text: str) -> List[str]:
    """The sequence each P line of a GFA 1.0 text spells: its segments (S, reverse-complemented on `-`) joined, each
    later one without its overlap with the previous one (the P line's overlap list, else the L line's CIGAR nM)."""
    seqs: Dict[str, str] = {}
    links: Dict[Tuple[str, str, str, str], Optional[int]] = {}
    paths: List[Tuple[List[str], List[str]]] = []
    for line in text.split("\n"):
        f = line.split("\t")
        if f[0] == "S" and len(f) >= 3:
            seqs[f[1]] = f[2]
        elif f[0] == "L" and len(f) >= 6:
            links[(f[1], f[2], f[3], f[4])] = _cigar_overlap(f[5])
        elif f[0] == "P" and len(f) >= 3:
            paths.append((f[2].split(","), f[3].split(",") if len(f) > 3 and f[3] != "*" else []))
    out = []
    for steps, ovl in paths:
        parts = []
        for j, st in enumerate(steps):
            name, o = st[:-1], st[-1]
            s = seqs.get(name, "")
            s = checks.revcomp(s) if o == "-" else s
            if j:
                k = _cigar_overlap(ovl[j - 1]) if j - 1 < len(ovl) else links.get(
                    (steps[j - 1][:-1], steps[j - 1][-1], name, o))
                s = s[k or 0:]
            parts.append(s)
        out.append("".join(parts))
    return out


def export_view(fmt: str, receipt: Optional[dict], text: Optional[str], what: Optional[str] = None) -> dict:
    """The comparable view of an export (FORMAT question Q3: the record holds the file): the receipt's records and
    the file's content in the shape W5 answers, so that `differences` compares a way's export with W5's."""
    out: dict = {}
    if isinstance(receipt, dict) and "records" in receipt:
        out["records"] = receipt["records"]
    if text is None:
        return out
    if fmt == "fasta":
        out["fasta"] = [{"id": r["id"], "sequence": r["sequence"]} for r in parse_fasta(text)]
    elif fmt == "gfa":
        ps = gfa_paths(text)
        out["gfa"] = {"n_paths": len(ps), "paths": sorted(ps)}
    elif fmt == "json" and what == "walk":
        doc = json.loads(text)
        out["arm"] = doc.get("arm")
        out["walks"] = doc.get("walks")
    elif fmt == "json" and what is None:
        out["json"] = normalise_full_result(json.loads(text))
    elif fmt in ("mgt", "save"):
        out["body_digest"] = body_digest(text)
        out["arms_digest"] = arms_digest(text)
    return out


def _walk_ids(d: _Doc, side: str, walks: Any) -> List[int]:
    n = len(d.arms[side]["paths"])
    ids = list(range(n)) if walks is None else ([walks] if isinstance(walks, int) else list(walks))
    if any(not isinstance(i, int) or isinstance(i, bool) or not 0 <= i < n for i in ids):
        raise _Refused("walks %r of %d" % (walks, n))
    return ids


def _graphlet_export(d: _Doc, args: dict, ctx: dict) -> dict:
    fmt = args.get("format", "fasta")
    what = args.get("what")
    if fmt not in ("fasta", "gfa", "json", "mgt") or (what is not None and fmt != "json"):
        raise _Refused("format %r / what %r" % (fmt, what))
    cons = ["/path", "/bytes", "/fields_cut", "/local", "/parsed", "/validation"]
    if fmt == "mgt":
        return _record({"body_digest": body_digest(d.mgt)},
                       {"/body_digest": "SPEC §7.5.2 the body (records without J and Z) of the document the tool "
                                        "holds: the standalone file carries it unchanged (FORMAT-v8 §3 body_digest)"},
                       cons + ["/records"], inputs=("mgt",))
    if fmt == "json" and what == "compare":
        return _declared("properties_only", "export json what=compare is graphlet_compare's lists")
    if fmt == "json" and what == "walk":
        if args.get("arm") is None or args.get("walks") is None:
            raise _Refused("what=walk without arm and walks", "bad_argument",
                           "DESIGN §6: what=walk with arm and walks writes their chains; the code is library-defined "
                           "(INVENTORY §2.13)")
        side = d.arm(args.get("arm"), required=True)
        ids = _walk_ids(d, side, args.get("walks"))
        walks, absent = [], {}
        for j, w in enumerate(ids):
            p = d.arms[side]["paths"][w]
            rows = _walk_rows(d, side, p, None)
            absent.update(_walk_row_extras(rows, "/walks/%d/segments" % j, None)[1])
            walks.append({"walk": w, "length_bp": p["length_bp"], "segments": rows})
        clauses = {"/arm": C_ARM_ECHO, "/walks/*/walk": C_PATHS, "/walks/*/length_bp": C_LENGTH}
        clauses.update({"/walks/*/segments/*/" + k: v.replace("first 8", "uncut").replace(
            "checked when at most 8 (INVENTORY §2.7)", "uncut") for k, v in _WALK_ROW_CLAUSES.items()})
        return _record({"arm": side, "walks": walks}, clauses,
                       cons + ["/records", "/walks/*/sequence", "/walks/*/segments/*/exact", "/handle", "/evidence"],
                       unordered=["/walks/*/segments/*/ends"], absent=absent)
    if fmt == "json":
        res = normalise_full_result(d.result)
        c_json = "SPEC §7.5.6: the full result is reproduced exactly from the text and its summary (the twin's " \
                 "result, normalised: DESIGN §9 T37, §2.5, SPEC §8.3)"
        timed = ["/json" + t for t in _time_measured(res)]
        for t in timed:   # a clock's value: left out of W5's answer, consistency-only (SPEC §6.8)
            parent, key = t.rsplit("/", 1)
            resolve_pointer({"json": res}, parent).pop(key.replace("~1", "/").replace("~0", "~"), None)
        clauses = {}
        for k in res:
            if k == "arms":
                for s, arm in res["arms"].items():
                    for ak in arm:
                        clauses["/json/arms/%s/%s" % (_esc(s), _esc(ak))] = c_json
            else:
                clauses["/json/" + _esc(k)] = c_json
        for t in timed:
            # a list or object holding a clock's value is compared field by field, that value left out
            entry = t.rsplit("/", 1)[0]
            container = entry.rsplit("/", 1)[0] if "/limitations/" in entry else entry
            if clauses.pop(container, None) is not None:
                clauses[container + ("/*/*" if container != entry else "/*")] = c_json
        return _record({"json": res}, clauses, cons + ["/records", "/json/timing", "/json/annotation/direct_reads",
                                                       "/json/annotation/keys_mapped"] + timed)
    sides = [d.arm(args.get("arm"), required=False)] if args.get("arm") is not None else d.sides
    if args.get("walks") is not None and args.get("arm") is None and len(d.sides) > 1:
        # which walks a list of ids names on two arms is defined by no clause (the library writes records for it)
        return _record({}, {}, cons + ["/records", "/fasta", "/gfa"], note="walks without an arm on two arms: "
                                                                        "library-defined, consistency-only")
    if fmt == "gfa":
        if args.get("arm") is not None or args.get("walks") is not None:
            return _record({}, {}, cons + ["/records", "/gfa"])
        mols = []
        for s in d.sides:
            fl = d.flanks(s)
            if fl is None:
                return _record({}, {}, cons + ["/records", "/gfa"])
            mols += [d.molecule(s, f) for f in fl]
        return _record({"gfa": {"n_paths": len(mols), "paths": sorted(mols)}},
                       {"/gfa/n_paths": "DESIGN §5 to_gfa: one P line per walk (both arms)",
                        "/gfa/paths": "DESIGN §9: the GFA paths spell the walks (with the seed, DESIGN §5 "
                                      "to_gfa(with_seed=True); " + C_SPELL + ")"},
                       cons + ["/records"])
    recs = []
    for s in sides:
        fl = d.flanks(s)
        if fl is None:
            return _record({}, {}, cons + ["/records", "/fasta"])
        for w in _walk_ids(d, s, args.get("walks")):
            recs.append({"id": "%s_%d" % (s, w), "sequence": d.molecule(s, fl[w])})
    return _record({"records": len(recs), "fasta": recs},
                   {"/records": "one FASTA record per walk exported (INVENTORY §2.13 records: the '>' count)",
                    "/fasta/*/id": "<side>_<path id> (INVENTORY §2.13, presentation), left arm first (SPEC §7.5.2 "
                                   "arm order), path id order",
                    "/fasta/*/sequence": "DESIGN §5 to_fasta(with_seed=True, orientation='natural'): the walk with the "
                                         "seed; " + C_SPELL},
                   cons + ["/fasta/*/header"])


def _graphlet_subtrie(d: _Doc, args: dict, ctx: dict) -> dict:
    sels = args.get("labels")
    sels = sels if isinstance(sels, list) else ([sels] if sels is not None else [])
    mode = args.get("mode", "any")
    if not sels or mode not in ("any", "all"):
        raise _Refused("labels %r / mode %r" % (sels, mode))
    side = d.arm(args.get("arm"), required=False) if args.get("arm") is not None else None
    ids = [d.resolve(s) for s in sels]
    view = {"selectors": [{"ref": d.labels[i]["ref"]} for i in ids], "mode": mode, "arm": side}
    sides = [side] if side else d.sides
    qualified = "for the selected labels"
    summ = {"view": dict(view), "labels": [d.label(i) for i in ids],
            "arms": {s: {"complete_to_bp": d.arms[s]["complete_to_bp"], "scope": d.arms[s]["completeness_scope"],
                         "qualified": qualified} for s in sides}}
    ev, ecl, econ = _evidence(d, sides)
    ev["view"] = dict(view)
    ev["qualified"] = qualified
    c_view = "DESIGN §5 subgraph(selectors, arm, mode): a view of the selected labels, selectors by ref (" + C_LABEL + ")"
    c_q = "DESIGN §5 subgraph: completeness is the backing complete_to_bp qualified 'for the selected labels'"
    clauses = {"/view/selectors": c_view, "/view/mode": c_view, "/view/arm": c_view,
               "/summary/view/selectors": c_view, "/summary/view/mode": c_view, "/summary/view/arm": c_view,
               "/summary/labels": c_view + ", in the order given",
               "/summary/arms/*/complete_to_bp": c_q + " (" + C_CTB + ")",
               "/summary/arms/*/scope": c_q + " (" + C_SCOPE + ")", "/summary/arms/*/qualified": c_q,
               "/evidence/view/selectors": c_view, "/evidence/view/mode": c_view, "/evidence/view/arm": c_view,
               "/evidence/qualified": c_q}
    clauses.update(ecl)
    econ = [x for x in econ if x not in ("/evidence/view", "/evidence/qualified")]
    return _record({"view": view, "summary": summ, "evidence": ev}, clauses,
                   econ + ["/handle", "/of", "/view/of", "/summary/view/of", "/evidence/view/of",
                           "/summary/arms/*/segments", "/summary/arms/*/walks", "/summary/outcome", "/fields_cut"])


def _graphlet_save(d: _Doc, args: dict, ctx: dict) -> dict:
    return _record({"body_digest": body_digest(d.mgt), "arms_digest": arms_digest(d.mgt)},
                   {"/body_digest": "SPEC §7.5.2: the saved file is H, J, the body unchanged (FORMAT-v8 §3 body_digest)",
                    "/arms_digest": "the saved file reads as the original (checks.graphlet_arms, SPEC §7.5)"},
                   ["/path", "/bytes", "/parsed", "/validation", "/fields_cut"], inputs=("mgt",))


def _graphlet_load(d: _Doc, args: dict, ctx: dict) -> dict:
    summ, clauses, cons = _summary_part(d, d.sides)
    cons = cons + ["/handle", "/fields_cut"]
    if ctx.get("max_bytes", 2048) < 16384:
        # the receipt drops its summary whole when it does not fit (DESIGN §6: the handle is always kept)
        return _record({}, {}, cons + ["/summary"])
    return _record({"summary": summ}, clauses, cons)


_OVERRIDE_FORBIDDEN = ("self", "g", "arm", "leaves", "reset_branches")


def _traverse_continue(d: _Doc, args: dict, ctx: dict) -> dict:
    """execute=False: the request a continuation of walk |walk| sends (SPEC §7.1 continuation; DESIGN §5.3). Checked:
    the seed (the continuation's bases; its labels in constrain mode, none in annotate mode), the request's loss
    budget lowered by the largest terminal loss of the continued labels, and its max_label_branches lowered by their
    largest terminal branch count ("unlimited" kept; the original with reset_branches). Whether the tool also states
    the loss_budget / branch_budget blocks is its presentation (INVENTORY §2.4: loss_budget?, branch_budget?), and how
    `overrides` merge into the request is library-defined: with overrides the strategy is consistency-only."""
    if _bool_arg(args, "execute", True):
        return _declared("not_applicable", "traverse_continue execute=true runs a backend request")
    side = d.arm(args.get("arm"), required=True)
    reset = _bool_arg(args, "reset_branches", False)
    ov = args.get("overrides")
    if ov is not None and (not isinstance(ov, dict) or set(ov) & set(_OVERRIDE_FORBIDDEN)):
        raise _Refused("overrides %r: an object that sets none of %s" % (ov, ", ".join(_OVERRIDE_FORBIDDEN)))
    w = args.get("walk")
    paths = d.arms[side]["paths"]
    if not isinstance(w, int) or isinstance(w, bool) or not 0 <= w < len(paths):
        raise _Refused("walk %r of %d" % (w, len(paths)))
    p = paths[w]
    c = p.get("continuation")
    if not c:
        raise _Refused("walk %d has no continuation" % w, None,
                       "SPEC §7.1: only a leaf ended by the radius, a resource reason or beam_pruned carries one; the "
                       "code is library-defined")
    strategy = ctx.get("strategy") or d.envelope.get("strategy") or {}
    cont_labels = sorted(c.get("labels") or [])
    terminal = {e["label"]: e for e in p.get("end_labels") or []}
    labels = strategy.get("labels") or {}
    branching = strategy.get("branching") or {}
    seed: dict = {"sequence": c.get("sequence") or ""}
    absent = {}
    if d.mode == "constrain":
        seed["labels"] = [d.labels[i]["name"] for i in cont_labels]
    else:
        absent["/request/seeds/0/labels"] = "SPEC §6.9, §5: in annotate mode seeds[].labels must be omitted (finding F3)"
    ans: dict = {"request": {"seeds": [seed]},
                 "parent": {"arm": side, "walk": w, "overlap_bp": len(c.get("sequence") or "")}}
    c_seed = "SPEC §7.1 continuation: valid /traverse input, a fresh seed whose labels are the continuation's " \
             "(names, SPEC §5 seeds[].labels; ascending label id)"
    clauses = {"/request/seeds/*/sequence": c_seed, "/request/seeds/*/labels": c_seed,
               "/parent/arm": C_ARM_ECHO, "/parent/walk": C_PATHS,
               "/parent/overlap_bp": "DESIGN §5.3: a new traversal with the overlap stated (the continuation's bases)"}
    cons = ["/request/strategy/*", "/request/seeds/*/seed_id", "/parent/handle", "/notes", "/ceiling_bytes",
            "/loss_budget", "/branch_budget", "/loss_budget_labels", "/fields_cut"]
    known = d.mode == "constrain" and all(i in terminal for i in cont_labels) and not ov
    if known:
        st: dict = {}
        lb = labels.get("loss_budget")
        if isinstance(lb, (int, float)) and not isinstance(lb, bool):
            largest = max((terminal[i].get("loss") or 0.0 for i in cont_labels), default=0.0)
            st["labels"] = {"loss_budget": lb - largest}
            clauses["/request/strategy/labels/loss_budget"] = \
                "SPEC §7.1 continuation: loss_budget lowered by the LARGEST terminal loss of the continued labels " \
                "(end_labels[].loss, the T records; DESIGN §5.3)"
        mlb = branching.get("max_label_branches")
        c_branch = "SPEC §7.1 continuation: max_label_branches lowered by the LARGEST terminal branch count of the " \
                   "continued labels, 'unlimited' kept, the original with reset_branches (DESIGN §5.3)"
        if mlb == "unlimited":
            st["branching"] = {"max_label_branches": "unlimited"}
        elif isinstance(mlb, int) and not isinstance(mlb, bool):
            largest_b = max((terminal[i].get("branches") or 0 for i in cont_labels), default=0)
            st["branching"] = {"max_label_branches": mlb if reset else max(mlb - largest_b, 0)}
        if "branching" in st:
            clauses["/request/strategy/branching/max_label_branches"] = c_branch
        if st:
            ans["request"]["strategy"] = st
    return _record(ans, clauses, cons, absent=absent,
                   note="overrides: their merge into the request is library-defined" if ov else None)


_TOOLS: Dict[str, Callable[[_Doc, dict, dict], dict]] = {
    "graphlet_summary": _graphlet_summary, "graphlet_walks": _graphlet_walks, "graphlet_walk": _graphlet_walk,
    "graphlet_labels": _graphlet_labels, "graphlet_splits": _graphlet_splits, "graphlet_sequence": _graphlet_sequence,
    "graphlet_export": _graphlet_export, "graphlet_subtrie": _graphlet_subtrie, "graphlet_save": _graphlet_save,
    "graphlet_load": _graphlet_load, "traverse_continue": _traverse_continue,
}


def answer(tool: str, args_resolved: dict, full_answer: dict, mgt_text: str, context: Optional[dict] = None) -> dict:
    """W5's record for one op (FORMAT-v8 §4); see the module docstring."""
    ctx = dict(context or {})
    args = {k: v for k, v in (args_resolved or {}).items() if k not in ("handle", "access_token", "cursor",
                                                                        "max_bytes", "budget")}
    if tool in PHASE_B or (tool == "graphlet_labels" and args.get("name") is not None):
        return _declared("not_implemented_phase_b", "%s: DESIGN-v8 §8 phase B" % tool)
    if tool in NOT_APPLICABLE:
        return _declared("not_applicable", "%s: a backend or store tool, outside graphlet processing" % tool)
    if tool == "graphlet_compare":
        return _declared("properties_only", "graphlet_compare: the compare_*_properties functions (DESIGN-v8 §2)")
    if tool not in _TOOLS:
        return _declared("not_applicable", "%s: not a graphlet tool" % tool)
    if ctx.get("view"):
        return _declared("view_consistency_only", "%s on a subtrie view: which walks and segments a view keeps is "
                                                  "not defined outside the library (DESIGN §5 subgraph)" % tool)
    twin = twin_difference(full_answer, mgt_text, ctx.get("seed_index", 0))
    if twin:
        return _declared("twin_mismatch", twin)
    d = _Doc(full_answer, mgt_text, ctx.get("seed_index", 0))
    try:
        return _TOOLS[tool](d, args, ctx)
    except _Refused as e:
        return _record(None, {}, ["/*"], declared="expected_error", expected_error=e.code,
                       note="%s (%s)" % (e, e.clause))


# ---------------------------------------------------------------------------------------------- compare properties
def _cut_walks(d: _Doc, side: str, depth: int) -> Optional[set]:
    """The arm's walks spelled in natural orientation and cut at |depth| outward bases (DESIGN §5.2.5: the DAG
    restricted to the comparison depth): right flank[:depth], left flank[-depth:]."""
    fl = d.flanks(side)
    if fl is None:
        return None
    return {f[:depth] if side == "right" else f[max(0, len(f) - depth):] for f in fl}


def _rows_of(ans: dict, side: str) -> List[dict]:
    return [r for r in ans.get("rows") or [] if r.get("side") == side]


def compare_self_properties(ans: dict, full_answer: Optional[dict], mgt_text: str) -> List[str]:
    """graphlet_compare(a, a): the verdict rules on a graphlet compared with itself (DESIGN §5 compare, §5.2.5,
    §3.1; SPEC §7.5.7): never `false` (one index identity, one oriented seed); without index_fp in H never a
    verdict that compares (true or qualified: "unverifiable, never equal"), and `unverifiable` only then; no
    difference under any verdict (counts zero or null, no rows); equal true exactly when comparable true; not true
    when the label evidence is a lower bound (DESIGN §5.2.5: what cannot be exact is qualified); neither true nor
    qualified for a side without bases in the modes that spell (claims, walks, prefix_subset: unknown). Where two
    rules give a non-comparing verdict (no bases and no index_fp) DESIGN states no precedence: either is accepted.
    -> problems."""
    d = _Doc(full_answer, mgt_text)
    probs = []
    comp = ans.get("comparable")
    if comp is False:
        probs.append("comparable false for a graphlet and itself (%s)" % ans.get("reason"))
    no_fp = d._h(14) is None
    if (no_fp and comp in (True, "qualified")) or (comp == "unverifiable" and not no_fp):
        probs.append("comparable %r with index_fp %s (SPEC §7.5.7, DESIGN §3.1)" % (comp, "absent" if no_fp else "set"))
    counts = ans.get("counts")
    if counts is not None and any(counts.get(k) for k in ("only_in_a", "only_in_b", "differ")):
        probs.append("differences %s between a graphlet and itself" % counts)
    if ans.get("rows"):
        probs.append("%d rows between a graphlet and itself" % len(ans["rows"]))
    if (ans.get("equal") is True) != (comp is True):
        probs.append("equal %r with comparable %r (DESIGN §5: equal only when comparable)" % (ans.get("equal"), comp))
    if comp is True and (d.result.get("outcome") or {}).get("label_evidence") not in (None, "complete"):
        probs.append("comparable true on label evidence %r (DESIGN §5.2.5: qualified)" %
                     d.result["outcome"]["label_evidence"])
    no_bases = any(d.flanks(s) is None for s in d.sides)
    if no_bases and ans.get("mode") in ("claims", "walks", "prefix_subset") and comp in (True, "qualified"):
        probs.append("comparable %r without bases (DESIGN §5.2.5: unknown)" % (comp,))
    return probs


def _swap_row(r: dict) -> dict:
    r = dict(r)
    r["side"] = {"a": "b", "b": "a"}.get(r.get("side"), r.get("side"))
    if r.get("side") == "differ" and "a" in r and "b" in r:
        r["a"], r["b"] = r["b"], r["a"]
    return r


def compare_swap_properties(ab: dict, ba: dict) -> List[str]:
    """graphlet_compare(a, b) vs (b, a), whole lists, modes claims, walks and labels: the same mode, verdict,
    equality and depth; only_in_a and only_in_b exchanged, differ kept; the rows the same up to the side exchange
    (DESIGN §5 compare: keyed by LabelRef and spelled walk, a symmetric relation reported per side). prefix_subset
    asks whether a is a prefix-subset of b (DESIGN §9 T39: tuned against exhaustive): directional, no swap
    property. -> problems."""
    if ab.get("mode") == "prefix_subset" or ba.get("mode") == "prefix_subset":
        return []
    probs = []
    for k in ("mode", "comparable", "equal", "depth_used"):
        if not same(ab.get(k), ba.get(k)):
            probs.append("%s: %r vs %r" % (k, ab.get(k), ba.get(k)))
    ca, cb = ab.get("counts"), ba.get("counts")
    if (ca is None) != (cb is None):
        probs.append("counts %r vs %r" % (ca, cb))
    elif ca is not None and (ca.get("only_in_a"), ca.get("only_in_b"), ca.get("differ")) != \
            (cb.get("only_in_b"), cb.get("only_in_a"), cb.get("differ")):
        probs.append("counts %r vs swapped %r" % (ca, cb))
    canon = lambda rows: sorted(json.dumps(r, sort_keys=True) for r in rows)   # noqa: E731
    if canon(_swap_row(r) for r in ab.get("rows") or []) != canon(ba.get("rows") or []):
        probs.append("rows differ beyond the side exchange (%d vs %d rows)" % (
            len(ab.get("rows") or []), len(ba.get("rows") or [])))
    return probs


def _compared_sides(da: _Doc, db: _Doc, arm: Any) -> List[str]:
    """The arms a compare covers: the one its `arm` argument names (the op's args_resolved: the compare answer itself
    names none), else every arm both retrievals walked."""
    want = _ARM_NAMES.get(arm) if isinstance(arm, str) else None
    return [s for s in ARMS if s in da.arms and s in db.arms and (arm is None or want == s)]


def compare_walks_properties(ans: dict, full_a: Optional[dict], mgt_a: str, full_b: Optional[dict], mgt_b: str,
                             arm: Optional[str] = None) -> List[str]:
    """mode walks, whole lists, a verdict that compares (true or qualified), on the arms the compare covers (|arm|:
    the op's `arm` argument): per arm, the walks only in a (b) are the spelled walks of a (b) cut at depth_used that
    b (a) lacks, and every differ row names a walk both have (DESIGN §5.2.5: the DAGs restricted to the comparison
    depth; walks spelled in natural orientation, SPEC §7.5.4). Rows name a walk by its bases (`walk`) and its arm.
    -> problems; [] when the verdict compares nothing."""
    if ans.get("comparable") not in (True, "qualified") or ans.get("mode") != "walks":
        return []
    D = ans.get("depth_used")
    if not isinstance(D, int) or isinstance(D, bool):
        return ["depth_used %r" % (D,)]
    da, db = _Doc(full_a, mgt_a), _Doc(full_b, mgt_b)
    probs = []
    sides = _compared_sides(da, db, arm)
    stray = [r for r in ans.get("rows") or [] if _ARM_NAMES.get(r.get("arm")) not in sides]
    if stray:
        probs.append("%d rows on an arm the compare does not cover (%s)" % (len(stray), sides))
    for side in sides:
        wa, wb = _cut_walks(da, side, D), _cut_walks(db, side, D)
        if wa is None or wb is None:
            probs.append("%s: no bases on a side under a comparing verdict" % side)
            continue
        got = {k: {r.get("walk") for r in _rows_of(ans, k) if _ARM_NAMES.get(r.get("arm")) == side}
               for k in ("a", "b", "differ")}
        if got["a"] != wa - wb:
            probs.append("%s: only_in_a %d walks, the twins' cut walks give %d" % (side, len(got["a"]), len(wa - wb)))
        if got["b"] != wb - wa:
            probs.append("%s: only_in_b %d walks, the twins' cut walks give %d" % (side, len(got["b"]), len(wb - wa)))
        if not got["differ"] <= (wa & wb):
            probs.append("%s: %d differ rows name a walk not in both" % (side, len(got["differ"] - (wa & wb))))
    return probs


def compare_prefix_subset_properties(ans: dict, full_small: Optional[dict], mgt_small: str,
                                     full_large: Optional[dict], mgt_large: str,
                                     arm: Optional[str] = None) -> List[str]:
    """mode prefix_subset on one seed at two radii, a = the smaller, on the arms the compare covers (|arm|: the op's
    `arm` argument): within the comparison depth the two retrievals are one trie (SPEC §6.10, T35: the walks equal
    the trie's at every depth <= complete_to_bp), so a comparing verdict has nothing to report (no rows, counts zero)
    and depth_used is at most the smaller retrieval's complete_to_bp on those arms (DESIGN §5.2.5). The twins' walks
    cut at depth_used must coincide on each covered arm: where they do not, the ENGINE broke SPEC §6.10's guarantee
    (a finding about the retrievals, named as such), and the property is not judged further. -> problems."""
    if ans.get("comparable") not in (True, "qualified") or ans.get("mode") != "prefix_subset":
        return []
    D = ans.get("depth_used")
    if not isinstance(D, int) or isinstance(D, bool):
        return ["depth_used %r" % (D,)]
    ds, dl = _Doc(full_small, mgt_small), _Doc(full_large, mgt_large)
    sides = _compared_sides(ds, dl, arm)
    probs = []
    ctb = min((ds.arms[s]["complete_to_bp"] for s in sides), default=0)
    if D > ctb:
        probs.append("depth_used %d beyond the smaller retrieval's complete_to_bp %d on %s" % (D, ctb, sides))
    for side in sides:
        ws, wl = _cut_walks(ds, side, D), _cut_walks(dl, side, D)
        if ws is None or wl is None or ws != wl:
            return probs + ["engine (SPEC §6.10, T35): the two radii's %s walks cut at %d differ, so they are no "
                            "single trie within the depth" % (side, D)]
    counts = ans.get("counts") or {}
    if any(counts.get(k) for k in ("only_in_a", "only_in_b", "differ")) or ans.get("rows"):
        probs.append("one trie within depth %d, yet counts %s and %d rows" % (D, counts, len(ans.get("rows") or [])))
    return probs
