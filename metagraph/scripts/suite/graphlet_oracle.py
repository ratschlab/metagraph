"""W5, the suite's independent oracle of the graphlet processing tools (manifest v8, DESIGN-v8 §2), phases A and B.

    answer(tool, args_resolved, full_answer, mgt_text, context=None) -> record
    compare(record, answer) -> {"outcome": agree | disagree | not_comparable | not_checked | inconclusive, "detail",
                                "checked", "differences"}
    is_cut_of(cut_value, full_value, cut_fields) / cut_problem(...)   # the suite's one cut rule
    python3 graphlet_oracle.py --run-dir RUN      # writes arms/oracle/graphlet_ops.jsonl (FORMAT v2 F4)

answer() gives the answer a library tool should give on the graphlet whose MGT text is |mgt_text| (the engine's
body, or W1's standalone document: answer() reads no J line, the run step only its strategy when it has no engine
body) and whose `detail: full` twin is |full_answer| (the direct-local response, or one result of it), for the
arguments W1 sent (FORMAT-v8 §2: W5 takes W1's args_resolved).
compare() is THE comparison of another way's answer with a W5 record (FORMAT v2 F2): the compare stage calls it.

The record (FORMAT-v8 §4 and FORMAT v2 F1):
    {"answer": {...} | None,
     "oracle": {"checked": {pointer: clause}, "consistency_only": [pattern], "inputs": [...],
                "unordered"?: [pointer], "position_ordered"?: [pointer], "absent"?: {pointer: clause},
                "rederived"?: [pointer], "declared"?: str, "expected_error"?: code | null, "note"?: str}}
- `answer` holds only what W5 derives; `checked` names each value by its RFC 6901 pointer with the SPEC/DESIGN
  clause it rests on (a clause saying "library-defined" or "INVENTORY §2 (presentation)" marks the few values whose
  only definition is the tool's own contract).
- Pointer patterns: in `consistency_only` (and accepted anywhere by compare()) a token `*` is any one list index or
  object key, e.g. "/rows/*/complete"; no tool answer has a key named `*`. Concrete pointers escape `~` and `/` as
  `~0` and `~1` (RFC 6901).
- `consistency_only`: the fields of the tool's answer W5 does not check (judged by W1..W4 agreeing only).
- `unordered`: checked lists compared as multisets (a tool list no clause orders). `position_ordered`: lists of
  events compared as sorted by position only (SPEC §7.1: "events at one at_bp are unordered"): non-decreasing at_bp,
  equal as multisets.
- `absent`: pointers (patterns expanded over the other answer) that must not exist there, each with its clause: an
  annotate continuation's seeds[].labels, merge_of or ends where SPEC says none, and every element past a whole
  list's end (the list's length, VO7).
- `rederived`: the checked pointers whose value W5 computes with a DESIGN formula from engine values, rather than
  reading an engine value (DESIGN-v8 §2: the report marks them "re-derived from DESIGN").
- `declared` says why a record checks no value: not_implemented_phase_b (no tool today: DESIGN-v8 §8's phase B
  is answered), properties_only (graphlet_compare: the compare_* functions below),
  not_applicable (backend and store tools), view_consistency_only (a tool on a subtrie view), cursor_page (a page
  resumed from an opaque cursor), path_unknown (a file path W5 cannot place inside or outside the export directory),
  local_limits (a budget under local limits): W5 makes no claim (compare: not_checked); twin_mismatch (the full twin is not the graphlet's walk: twin_difference, DESIGN-v8 §2's
  twin pre-check), no_input (the run step found no twin or document), unresolved_handle (a handle symbol naming no
  scenario or op), selector_error (W1's arguments did not resolve), unreadable_w1 (a torn W1 line): W5 could not
  judge (compare: inconclusive, which blocks); and expected_error: the arguments name no answer, or the full twin
  states the seed failed (SPEC §7.1: no graphlet; every op that reads one expects no_graphlet), the tool must
  refuse; then
  `expected_error` is the code W5 expects (null: any code), the note gives why and the clause (most codes are the
  library's vocabulary, so "library-defined").
- Cut rows: a field the other answer marks cut to fit its page (FORMAT v2 A3: `_cut_fields` on each row of a page
  with row_truncated) must be a cut of W5's value (REVIEW D-R8: a prefix; a sub-multiset of a set).
- Fitting errors (FIT_ERRORS: result_too_large and kin) are not_comparable against W5 (FORMAT v2 F3).

`context` (all optional): seed_index (0); paging ("all", the default: lists are the whole list, every page
concatenated; otherwise rows are checked up to `rows`, the reference way's row count, and are consistency-only
without it); max_bytes (2048) and sequence_max_bytes (16384), the configured ceilings (DESIGN-v8 §3: from the
capabilities' limits block); view (the op runs on a subtrie handle); strategy (the graphlet envelope's normalized
strategy, for traverse_continue and the twin's memory budget below; default the full twin's).

Independence (DESIGN-v8 §2, REVIEW-DESIGN-v8 D-R5): this module imports the standard library and the suite's own
`checks` (its MGT reader, stdlib only) and nothing of `metagraph`; it cites SPEC-core
(docs/SPEC-labeled-traversal-core.md, "SPEC") and DESIGN (docs/DESIGN-traverse-graphlet.md) only. Every value comes,
by preference, from the engine's `full` twin (segments, paths with end_labels and continuation, runs, splits,
label_summary, limitations: computed by the C++ walker) and the MGT tokens (H, S); a field whose only definition is
a tool's own presentation is cited "INVENTORY §2 (presentation)" (list caps, tie orders, string forms), and a field
with no definition outside the library, or only DESIGN pseudocode the twin does not compute, is consistency-only.

Values a clock measured (a time limitation's observed, a time resource stop's amounts) differ between the twin's
execution and the graphlet's (SPEC §6.8) and are never checked.

A twin without the graphlet request's memory budget (VERIFY-D VD4: v8 sends the full twin without
bounds.max_memory_mb, since `full` detail under a memory budget keeps a merge's arrival order, SPEC §7.1 R21 (4); the
budget is per seed, §6.8): the graphlet answer states one seed-level memory_bound_soft the twin does not (SPEC §7.0:
stated exactly when the knob is set). W5 then reads the seed-level entries of graphlet_summary's limitations and of
the json export from the graphlet's own K records (SPEC §7.5.2 K: the JSON entries, seed level first, in order;
VERIFY-O-VD4 suggestion 1), checks everything else as before, and names the difference in every record's
oracle.inputs (_memory_gap: only when the twin's strategy has no budget and its result states no memory_bound_soft).

Coverage (phases A and B; [S] SPEC, [D] DESIGN, [P] presentation per INVENTORY §2, [L] library-defined (an expected
error's code); [R] re-derived from DESIGN (oracle.rederived); c/o = consistency-only). Every whole list's length is
checked (VO7); refusals are expected errors with their code [L].
  every answer's evidence: complete_to_bp, scope [S §7.4], support, reconverge [S §7.5.2 H], outcome [S §7.0],
      exact false where DESIGN §6's condition fails (label evidence not complete, a cut list) | c/o: limitations
      (the tool's digest, finding F1), exact otherwise, coordinates, view.of
  graphlet_summary summary: seed.length_bp [S §7.1], seed.mode [S §7.4 label_mode], outcome [S §7.0],
      index.ns, index.verifiable [S §7.5.2 H, §7.5.7], arms.*: status, complete_to_bp, scope [S §7.4], walks,
      segments, splits, merges, bases, max_bp [S §7.5.5 counts, §7.5.2 A] | c/o: seed.seed_id, seed.labels,
      top_labels, caveats, cut flags
  graphlet_summary limitations: outcome, resource_stop [S §7.0], rows (each entry, arm) and total [S §7.0, §7.5.2 K
      order]; under a twin without the memory budget the seed-level rows are the graphlet's K records [S §7.5.2 K]
      (VD4) | c/o: rows.*.informational; a clock's values (SPEC §6.8); rows with arm= beside seed-level entries
  graphlet_walks: arm, total, rows.*.walk (in the rank's order: support [D §5][R], length, id [P ties]), length_bp,
      end [S §7.5.6, §7.1; string P], n_alive [S §7.5.6], n_full, labels_full (+_truncated) [D §5.1 via
      end_labels.route_bp and runs.from_bp][R], sequence (spell tail/full) [S §7.5.4], label filter (alive at the
      leaf) and route_consistent (filtered.merge_entered: entered through a non-first parent) [D §6, §5.1][R];
      coordinates=true without record coordinates: expected no_coordinates [S §7.1][L] | c/o: complete,
      beyond_certified_bp, coordinates, hint, filtered.not_from_seed; rank=loss order (no key in DESIGN);
      annotate: n_alive, n_full, labels_full, rank=support order, label filters
  graphlet_walk: arm, walk, length_bp, total, rows.*: segment, from_bp, to_bp, labels_in(_total), labels_out(_total)
      [S §7.1, §7.4], in, out [S §7.1; first 8 P], merge_of (absent on one parent) [S §7.5.6], ends (a set while
      ENDS_UNORDERED, G3; absent without one) [S §7.5.6 label_end events], continuation {length_bp, loss_used,
      labels} (absent without one) [S §7.1] | c/o: rows.*.exact; ends beyond 8; labels_out_total on cut lists
  graphlet_labels (table): total, rows (name, ref, arm, direct_bp, reach_bp, reentries, runs) [S §7.1
      label_summary], order (-rank, ref, arm) [D §6; ties P], min_direct_bp, at_bp (constrain: a run holds the base,
      S §7.1; annotate: a recorded label set holds it, S §7.4) | c/o: at_bp where a cut label set covers the base
  graphlet_splits: total, rows (at_bp, segment, kind, labels_before, branches.char, branches.labels =
      labels_distinct) in (at_bp, first child) order [S §7.4, §7.5.6], min_labels_before
  graphlet_sequence: arm, orientation, with_seed, from, length, sequence [S §7.5.4], to (an uncut slice) | c/o:
      ceiling_bytes, truncated/next_from/to of a cut slice, left-arm walk orientation with the seed
  graphlet_export: fasta records, ids and sequences [D §5 to_fasta, S §7.5.4]; gfa P count and the walks the paths
      spell [D §9]; json = the normalised full result per key [S §7.5.6, D §9 T37] (seed-level limitations: the
      graphlet's K records under a twin without the memory budget [S §7.5.2 K], VD4); json what=walk chains (arm and
      walks required [L §2.13]); mgt body digest [S §7.5.2, FORMAT §3]; fasta without bases: an expected error [L] | c/o:
      header tags, gfa S/L lines, gfa with arm or without bases, records of gfa/json/mgt, what=compare, fasta walks
      without arm
  graphlet_subtrie: view selectors, mode, arm, summary labels, per-arm complete_to_bp, scope, qualified [D §5
      subgraph] | c/o: the view's walk and segment counts; every tool on the view handle
  graphlet_save / graphlet_load: the saved text's body digest and its graphlet_arms (round trip) [S §7.5.2];
      load's summary as graphlet_summary's
  traverse_continue (execute=False): request.seeds (sequence; labels, absent in annotate mode [S §6.9, §5;
      finding F3]), the request's loss_budget and max_label_branches ('unlimited' kept, the original with
      reset_branches), parent arm/walk/overlap_bp [S §7.1 continuation, D §5.3] | c/o: the loss_budget /
      branch_budget blocks [P §2.4], labels.extra, notes, the rest of the strategy, the whole strategy with
      overrides (their merge [L]); execute=True (backend)
  graphlet_compare: properties only (compare_self_properties, compare_swap_properties, compare_walks_properties,
      compare_prefix_subset_properties; each takes the compare's arm) [D §5, §5.2.5, §3.1, §9 T39, S §6.10, §7.5.7]
  graphlet_claims (phase B): one claim per run end in run id order [S §6.9]: label, from_bp, to_bp, entered_by,
      from_label, cost, loss, branches, reason [S §7.1 runs[], §7.5.2 R], qualifier [S §7.5.1], path_id of a run
      ending at a leaf [S §7.5.6], kind boundary of a zero-length run [D §5.1], route_bp where R's stamp is the
      displayed one [S §7.1; PB2], evidence_from [D §5.1][R] (null for a boundary claim, PB1), total and filtered
      (labels [D §5], route_consistent: merge_entered where the SPEC §7.1 stamp is above 0, T at a leaf, R inside
      [D §6, S §7.1; F2, VO-B2], min_bp [P §2.11][R]) | c/o: other kinds, end_class, support, exact, path_id of
      inner runs, route_bp where R and T differ (PB2), the selection where R's and the displayed stamp disagree,
      filtered with min_bp > 0 (PB3), labels=[], coordinate columns, hint, caveat; annotate claims (D §5.2.1)
  graphlet_labels name= (phase B): the label's runs per arm (left first) in run id order: arm, from_bp, to_bp [S
      §7.1], end (the qualifier's text or the enum [S §7.5.6]), evidence_from [D §5.1][R] (the formula also for a
      zero-length run, PB4), merged_into [S §7.5.2 R m, D §4], route [D §5 routes][R]; name, ref, total; unknown /
      ambiguous name [L] | c/o: walks_below (undefined, VO-B3), the end of Lw and m runs, coordinate columns;
      annotate (D §5.2.1 witness route)
  graphlet_support (phase B): total, rows from_bp, to_bp, n_labels, n_labels_total, joined, left after the first row:
      constrain, each run's displayed evidence from its anchor, kept on the walk [D §5, §5.1][R] (G4: the
      library's presence profile is a defect);
      annotate, the recorded label sets along the walk [S §7.4] | c/o: the first row's lists, why, exact, step=, who
      joined or left across a cut list, n_labels_total on a constrain arm with cut lists
"""

from __future__ import annotations

import argparse
import copy
import gzip
import hashlib
import json
import os
import posixpath
import re
import sys
import time
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import checks
import suite_common

ARMS = ("left", "right")
_ARM_NAMES = {"left": "left", "right": "right", "l": "left", "r": "right"}
_SUPPORT = {"k": "kmer", "t": "trace"}            # SPEC §7.5.2 H <support k|t>
_RECONVERGE = {"m": "merge", "k": "keep"}         # SPEC §7.5.2 H <reconverge m|k>
# G3 (VO9, F5): SPEC §7.1 / DESIGN §2.5 order a walk row's `ends` by position; the library does since gr/v8-libfix
# (d39b6b38b), so W5 compares them position-ordered (True would compare them as a set, the pre-fix reading)
ENDS_UNORDERED = False
LIST_CAP = 8      # the label lists a row shows (labels_full, in, out, ends): INVENTORY §2.6/§2.7 (presentation)
# bytes of a sequence page beyond the fields W5 answers: the handle (a UUID at most), ceiling_bytes, truncated and
# next_from with their keys, and the JSON frame; a slice fits when it, W5's own page without the slice and this
# allowance stay within the sequence ceiling (graphlet_sequence)
_SEQ_PAGE_EXTRA = 512

PHASE_B: Tuple[str, ...] = ()   # every phase-B tool is answered (graphlet_support, graphlet_claims, labels name=)
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
C_NO_BASES_EXPORT = "library-defined: the tool refuses a FASTA of a retrieval whose strings SPEC §7.1 omits " \
                    "(sequences: false), with a code of its own"
C_WHOLE_LIST = "the whole list W5 derives (every page assembled, FORMAT-v8 §3): no element beyond it"
C_BAD_ARG = "library-defined: an argument outside the tool's domain is bad_argument (INVENTORY §1.1)"
C_MEMORY_TWIN = "SPEC §7.1 (R21 (4): `full` detail under a memory budget keeps a merge's arrival order, so the full " \
                "twin runs without bounds.max_memory_mb) and §6.8 (max_memory_mb is a per-seed budget); SPEC §7.0: " \
                "memory_bound_soft is one seed-level entry, stated exactly when bounds.max_memory_mb is set"


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
        self.mode = self.result.get("label_mode") or {"c": "constrain", "a": "annotate"}.get(self.h(6), self.h(6))
        self.labels = self._label_list()

    def h(self, i: int) -> Optional[str]:
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
        """Per path id, the flank in natural orientation, or None without sequences (SPEC §7.1: `sequence` holds the
        bases a segment adds, natural orientation; a right flank is the concatenation root -> leaf, a left flank
        leaf -> root)."""
        a = self.arms[side]
        segs = {sg["id"]: sg for sg in a.get("segments") or []}
        out = []
        for p in a.get("paths") or []:
            parts = [segs.get(i, {}).get("sequence") for i in p.get("segments") or []]
            if any(x is None for x in parts):
                return None
            out.append("".join(parts) if side == "right" else "".join(reversed(parts)))
        return out

    def molecule(self, side: str, flank: str) -> str:
        """The walk with the seed, natural orientation: natural(left) + seed or seed + natural(right) (SPEC §7.5.4)."""
        return flank + self.seed if side == "left" else self.seed + flank

    # ------------------------------------------------------------------ lineages (phase B)
    def r_records(self, side: str) -> List[List[str]]:
        """The arm's R records in order (SPEC §7.5.2: run id = ordinal, aligned with the twin's runs[]): R <segment>
        <label> <from_bp> <to_bp> <end> <route_bp> <from_label:cost|*> <prev_run|*> <structural_successors|*>
        <branches> <loss> [<needed_budget>]."""
        if not hasattr(self, "_r"):
            self._r: Dict[str, List[List[str]]] = {"left": [], "right": []}
            cur = None
            for line in (self.mgt or "").split("\n"):
                f = line.split(" ")
                if f[0] == "A" and len(f) > 1:
                    cur = "left" if f[1] == "l" else "right"
                elif f[0] == "R" and cur is not None:
                    self._r[cur].append(f)
        return self._r.get(side, [])

    def runs(self, side: str) -> List[dict]:
        return self.arms[side].get("runs") or []

    def aligned(self, side: str) -> bool:
        """Whether the R records and the twin's runs[] are one run each, in one order (label, from_bp, to_bp)."""
        r, runs = self.r_records(side), self.runs(side)
        return len(r) == len(runs) and all(len(f) > 8 and f[2] == str(x["label"]) and f[3] == str(x["from_bp"])
                                           and f[4] == str(x["to_bp"]) for f, x in zip(r, runs))

    def prev_run(self, side: str, i: int) -> Optional[int]:
        """The run a switch-entered run continues (SPEC §7.5.2 R prev_run; * = none)."""
        f = self.r_records(side)[i]
        return None if f[8] == "*" else int(f[8])

    def end_code(self, side: str, i: int) -> str:
        """The run's end token (SPEC §7.5.2 R end: code[+qualifier] | Lw | m)."""
        return self.r_records(side)[i][5]

    def segment(self, side: str, sid: int) -> dict:
        return self.arms[side]["segments"][sid]

    def chain(self, side: str, sid: int) -> List[int]:
        """The first-parent chain root -> segment |sid| (SPEC §7.5.6: the displayed chain through parents[0])."""
        out = [sid]
        while self.segment(side, out[-1]).get("parents"):
            out.append(self.segment(side, out[-1])["parents"][0])
        return out[::-1]

    def label_at(self, side: str, i: int, depth: int) -> int:
        """The INCOMING label of run |i|'s lineage at outward depth |depth|, before any switch committed there
        (DESIGN §5.1: walk R.prev_run back across switches while the run begins at or after the depth)."""
        runs = self.runs(side)
        while runs[i]["from_bp"] >= depth and self.prev_run(side, i) is not None:
            i = self.prev_run(side, i)
        return runs[i]["label"]

    def evidence_from(self, side: str, i: int, merge_part: bool = False) -> Optional[int]:
        """DESIGN §5.1 displayed-path evidence of run |i| (pseudocode, re-derived): p = the first-parent chain from the
        root to the run's anchor segment, t = to_bp; the latest merge s on p with from_bp(s) <= t whose first
        partition (the twin's labels_via_parent[0], SPEC §7.1) lacks the lineage's incoming label at from_bp(s); at
        least the run's own from_bp. None where a merge on p states no partitions."""
        run = self.runs(side)[i]
        ef = 0
        merges = [self.segment(side, sid) for sid in self.chain(side, run["segment"])
                  if len(self.segment(side, sid).get("parents") or []) > 1]
        for sg in sorted(merges, key=lambda x: -x["from_bp"]):
            if sg["from_bp"] > run["to_bp"]:
                continue
            parts = sg.get("labels_via_parent")
            if not isinstance(parts, list) or not parts:
                return None
            if self.label_at(side, i, sg["from_bp"]) not in parts[0]:
                ef = sg["from_bp"]
                break
        return ef if merge_part else max(ef, run["from_bp"])

    def route(self, side: str, i: int) -> Optional[List[int]]:
        """Run |i|'s own route, root -> its anchor segment (DESIGN §5 routes: through merges via the partitions,
        re-derived): at a merge, the parent whose partition (labels_via_parent, SPEC §7.1) holds the lineage's
        incoming label there. None where a merge states no partition holding it."""
        out = [self.runs(side)[i]["segment"]]
        while True:
            sg = self.segment(side, out[-1])
            parents = sg.get("parents") or []
            if not parents:
                return out[::-1]
            if len(parents) == 1:
                out.append(parents[0])
                continue
            parts = sg.get("labels_via_parent") or []
            lab = self.label_at(side, i, sg["from_bp"])
            k = next((k for k, part in enumerate(parts) if lab in part), None)
            if k is None or k >= len(parents):
                return None
            out.append(parents[k])


# ---------------------------------------------------------------------------------------------- the record
def _esc(key: Any) -> str:
    return str(key).replace("~", "~0").replace("/", "~1")


def _matches(pattern: str, pointer: str) -> bool:
    p, q = pattern.split("/"), pointer.split("/")
    return len(p) == len(q) and all(a == "*" or a == b for a, b in zip(p, q))


def _record(ans: Optional[dict], clauses: Dict[str, str], consistency: Sequence[str],
            inputs: Sequence[str] = ("full", "mgt"), unordered: Sequence[str] = (), absent: Any = (),
            position_ordered: Sequence[str] = (), partial: Sequence[str] = (), rederived: Sequence[str] = (),
            **extra: Any) -> dict:
    """The W5 record: every value of |ans| claimed by the first clause pattern that matches its pointer (a pattern
    naming an object or a list claims it whole). A value no clause claims is a bug of this module: W5 answers only
    what it checks (FORMAT-v8 §4). |unordered|: patterns of checked lists compared as multisets (a tool list whose
    order no clause fixes); |absent|: pointers (with their clause) that must not exist in the other answer, e.g. a
    seed's labels in annotate mode (FORMAT questions Q2, Q4). A list W5 claims element by element is a whole list
    unless named in |partial| (rows of a listing W1 did not read to its end) or consistency-only: its length is
    checked by an absent pointer at its end (VO7: an extra or repeated element is a difference)."""
    checked: Dict[str, str] = {}
    lengths: Dict[str, str] = {}

    def visit(node: Any, ptr: str) -> None:
        for pat, clause in clauses.items():
            if _matches(pat, ptr):
                checked[ptr] = clause
                return
        if isinstance(node, list):
            end = "%s/%d" % (ptr, len(node))
            if ptr not in partial and not any(_matches(c, end) for c in consistency):
                lengths[end] = C_WHOLE_LIST
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
    posord = sorted(p for p in checked if any(_matches(u, p) for u in position_ordered))
    if posord:
        oracle["position_ordered"] = posord
    red = sorted(p for p in checked if any(_matches(u, p) for u in rederived))
    if red:
        oracle["rederived"] = red
    absent = dict(absent, **lengths)
    if absent:
        oracle["absent"] = absent
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


# JSON value equality as the suite compares parsed values (booleans only with booleans, numbers by value)
same = suite_common.json_equal


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


# The compare stage's error codes of fitting (an answer that did not fit a ceiling or a service lane, not a statement
# about the graphlet): `not_comparable` against W5, which makes no claim about fitting (FORMAT v2 F3, VO5)
FIT_ERRORS = ("result_too_large", "receipt_too_large", "too_large", "body_too_large", "op_busy", "op_timeout",
              "op_memory", "local_budget_exceeded")
# declarations by which W5 makes no claim about the call, by design: compare() answers not_checked
NO_CLAIM = ("not_implemented_phase_b", "properties_only", "not_applicable", "view_consistency_only", "cursor_page",
            "path_unknown", "local_limits")
# declarations by which W5 could not judge a call it should judge (an input missing or inconsistent, a handle or a
# W1 line it cannot read): compare() answers inconclusive, which blocks (R2-5; owner policy 2026-10-14)
CANNOT_JUDGE = ("no_input", "unresolved_handle", "twin_mismatch", "selector_error", "unreadable_w1")


def expand(pattern: str, doc: Any) -> List[str]:
    """The concrete pointers of |doc| that a pointer pattern names: a token `*` is any list index or object key (no
    tool answer has a key named `*`); a pattern without `*` is itself, whether or not it resolves."""
    if "*" not in pattern.split("/"):
        return [pattern]
    out: List[str] = []

    def walk(node: Any, toks: List[str], ptr: str) -> None:
        if not toks:
            out.append(ptr)
            return
        t, rest = toks[0], toks[1:]
        if t == "*":
            items = node.items() if isinstance(node, dict) else enumerate(node) if isinstance(node, list) else ()
            for k, v in items:
                walk(v, rest, ptr + "/" + _esc(k))
            return
        key = t.replace("~1", "/").replace("~0", "~")
        if isinstance(node, dict) and key in node:
            walk(node[key], rest, ptr + "/" + t)
        elif isinstance(node, list) and t.isdigit() and int(t) < len(node):
            walk(node[int(t)], rest, ptr + "/" + t)

    walk(doc, pattern.split("/")[1:], "")
    return out


def _cut_by_other(other: Any, ptr: str) -> bool:
    """Whether |other| names the field at |ptr| as cut to fit its page (DESIGN §6: a row larger than max_bytes comes
    alone with row_truncated and the fields cut named): FORMAT v2 A3 marks each row of such a page with
    `_cut_fields`, read on the object that holds the field (a row's own row_truncated + cut_fields too)."""
    parent, key = ptr.rsplit("/", 1)
    try:
        holder = resolve_pointer(other, parent)
    except KeyError:
        return False
    if not isinstance(holder, dict):
        return False
    return key in (holder.get("_cut_fields") or ()) or bool(holder.get("row_truncated")
                                                            and key in (holder.get("cut_fields") or ()))


def _sub_multiset(small: list, big: list) -> bool:
    pool = [_canon(x) for x in big]
    for x in small:
        if _canon(x) not in pool:
            return False
        pool.remove(_canon(x))
    return True


def _at(x: Any) -> Any:
    return x.get("at_bp", 0) if isinstance(x, dict) else 0


def _position_equal(mine: list, theirs: list) -> bool:
    """A list of events sorted by position only (SPEC §7.1: "events at one at_bp are unordered"): |theirs| in
    non-decreasing at_bp and, as a multiset, |mine|."""
    return all(_at(a) <= _at(b) for a, b in zip(theirs, theirs[1:])) and len(mine) == len(theirs) \
        and _sub_multiset(theirs, mine)


def _is_cut_of(mine: Any, theirs: Any, order: str) -> bool:
    """|theirs| a cut of |mine| in the suite's words (REVIEW D-R8): a string or a list a prefix; an unordered list a
    sub-multiset; a position-ordered list a sub-multiset in position order that keeps the smallest positions."""
    if isinstance(mine, str) and isinstance(theirs, str):
        return mine.startswith(theirs)
    if isinstance(mine, list) and isinstance(theirs, list):
        if order == "unordered":
            return _sub_multiset(theirs, mine)
        if order == "position":
            ats = sorted(_at(x) for x in mine)[:len(theirs)]
            return _sub_multiset(theirs, mine) and [_at(x) for x in theirs] == ats
        return len(theirs) <= len(mine) and all(same(x, y) for x, y in zip(mine, theirs))
    return same(mine, theirs)


_CUT_MARKERS = ("_cut_fields", "row_truncated", "cut_fields")


def cut_problem(cut_value: Any, full_value: Any, cut_fields: Sequence[str] = (), unordered: Sequence[str] = (),
                position_ordered: Sequence[str] = ()) -> Optional[str]:
    """The suite's one cut rule (REVIEW D-R8, DESIGN §6 row_truncated / cut_fields; FORMAT v2 A3): why |cut_value|
    is NOT a cut of |full_value|, or None when it is. Used by differences() and by the compare stage's cut-row checks
    across ways (lane C).

    - Two objects (a row, or a page base): every key of |full_value| is in |cut_value| (the cut markers
      `_cut_fields`, `row_truncated`, `cut_fields` are ignored on both sides); a key named in |cut_fields| holds a cut
      of the full value, every other key the same value (json_equal).
    - A string: a prefix of the full one. A list: a prefix; for a field named in |unordered|, a sub-multiset; in
      |position_ordered| (events sorted by position only, SPEC §7.1), a sub-multiset in non-decreasing at_bp that
      keeps the smallest positions. A scalar is never cut (INVENTORY §1.2): equal. A cut may be empty (a string
      cut to "", a list to []): the tool's row cut shortens a field as far as the page needs, to nothing if need be.
    |unordered| / |position_ordered| name fields of an object, or, for a bare list, use "" for the list itself."""
    def order(key: str) -> str:
        return "unordered" if key in unordered else "position" if key in position_ordered else "ordered"

    if isinstance(full_value, dict) and isinstance(cut_value, dict):
        for k, v in full_value.items():
            if k in _CUT_MARKERS:
                continue
            if k not in cut_value:
                return "%s: missing in the cut" % k
            if k in cut_fields:
                if not _is_cut_of(v, cut_value[k], order(k)):
                    return "%s: not a cut of the full value (%s)" % (k, order(k))
            elif not same(v, cut_value[k]):
                return "%s: differs and is not named cut" % k
        extra = [k for k in cut_value if k not in full_value and k not in _CUT_MARKERS]
        return "%s: only in the cut" % extra[0] if extra else None
    if not _is_cut_of(full_value, cut_value, order("")):
        return "not a cut of the full value (%s)" % order("")
    return None


def is_cut_of(cut_value: Any, full_value: Any, cut_fields: Sequence[str] = (), unordered: Sequence[str] = (),
              position_ordered: Sequence[str] = ()) -> bool:
    """Whether |cut_value| is a cut of |full_value| (cut_problem() is None): the suite's one cut rule."""
    return cut_problem(cut_value, full_value, cut_fields, unordered, position_ordered) is None


def differences(record: dict, other: Optional[dict]) -> List[dict]:
    """Where |other| (another way's assembled answer, service fields removed) departs from W5: [{pointer, clause,
    oracle, other | missing | present}]. Every checked pointer (a pattern expanded over W5's answer) must resolve in
    |other| with the same value; lists named `unordered` compare as multisets, `position_ordered` ones as events
    sorted by position only; a field |other| marks cut to fit (FORMAT v2 A3 `_cut_fields`) must be a cut of W5's
    (_is_cut_of); a pointer named `absent` (a pattern expanded over |other|) must not resolve. A record declaring
    `expected_error` requires an error in |other|, of the code W5 states (oracle.expected_error, when not null);
    any other declaration checks nothing. compare() is the entry the compare stage calls."""
    o = record["oracle"]
    if o.get("declared") == "expected_error":
        code = error_code(other)
        if code is None:
            return [{"pointer": "/error", "clause": "DESIGN §6: errors are structured results (%s)" % o.get("note"),
                     "oracle": o.get("expected_error") or "an error", "missing": True}]
        if o.get("expected_error") and code != o["expected_error"]:
            return [{"pointer": "/error/code", "clause": o.get("note"), "oracle": o["expected_error"], "other": code}]
        return []
    ans = record.get("answer")
    unord = {q for p in o.get("unordered", ()) for q in expand(p, ans)}
    posord = {q for p in o.get("position_ordered", ()) for q in expand(p, ans)}
    out = []
    for pat, clause in o["checked"].items():
        for ptr in expand(pat, ans):
            mine = resolve_pointer(ans, ptr)
            try:
                theirs = resolve_pointer(other, ptr)
            except KeyError:
                out.append({"pointer": ptr, "clause": clause, "oracle": mine, "missing": True})
                continue
            order = "unordered" if ptr in unord else "position" if ptr in posord else "ordered"
            if _cut_by_other(other, ptr):
                if not _is_cut_of(mine, theirs, order):
                    out.append({"pointer": ptr, "clause": clause + "; cut to fit (DESIGN §6, REVIEW D-R8): a prefix",
                                "oracle": mine, "other": theirs})
                continue
            if isinstance(mine, list) and isinstance(theirs, list) and order != "ordered":
                ok = _position_equal(mine, theirs) if order == "position" else \
                    same(sorted(mine, key=_canon), sorted(theirs, key=_canon))
            else:
                ok = same(mine, theirs)
            if not ok:
                out.append({"pointer": ptr, "clause": clause, "oracle": mine, "other": theirs})
    for pat, clause in (o.get("absent") or {}).items():
        for ptr in expand(pat, other):
            try:
                theirs = resolve_pointer(other, ptr)
            except KeyError:
                continue
            out.append({"pointer": ptr, "clause": clause, "oracle": "absent", "present": theirs})
    return out


def compare(record_w5: dict, answer: Any) -> dict:
    """THE comparison of one way's answer with a W5 record (FORMAT v2 F2; lane C's compare_w5 calls it).

    |record_w5|: answer()'s record or a W5 line of arms/oracle/graphlet_ops.jsonl (both carry `answer` and `oracle`).
    |answer|: the other way's per-call record (FORMAT §3: its `answer`, or its `error` {"code", "fields"}, B2) or a
    bare tool answer (an error is {"error": "<code>", ...}), service fields already removed (DESIGN-v8 §5). For
    graphlet_export and graphlet_save the answer to pass is export_view(format, receipt, the file's text) (FORMAT v2
    B3: the file kept under the way's raw/), the shape W5 answers for a file.
    -> {"outcome": "agree" | "disagree" | "not_comparable" | "not_checked" | "inconclusive", "detail": str,
        "checked": n, "differences": [...]}:
    - not_checked: W5 makes no claim about the call by design (NO_CLAIM: phase B, properties only, not applicable,
      a view; or a record that checks nothing: consistency-only);
    - inconclusive (blocks): W5 could not judge a call it should (CANNOT_JUDGE: no input, an unresolved handle
      symbol, a twin mismatch, a W1 selector error, an unreadable W1 line), or the record is no W5 record;
    - not_comparable: the other way answered a fitting error (FIT_ERRORS: F3), about which W5 says nothing;
    - otherwise differences() decides: an expected error must be one (of W5's code when stated); an answer W5
      derives must hold every checked value (cut fields as cuts, unordered and position-ordered lists as such) and
      none of the absent keys; an error where W5 answers disagrees."""
    o = record_w5.get("oracle") or {}
    if isinstance(answer, dict) and "answer" in answer and ("error" in answer or "args_resolved" in answer):
        other = {"error": answer["error"]} if answer.get("error") else answer.get("answer")
    else:
        other = answer
    declared = o.get("declared")
    n = len(o.get("checked") or {})
    if not isinstance(record_w5.get("oracle"), dict) or not isinstance(o.get("checked"), dict) or \
            (declared and declared not in NO_CLAIM + CANNOT_JUDGE + ("expected_error",)):
        return {"outcome": "inconclusive", "detail": "not a W5 record W5 can be judged by (declared %r)" % declared,
                "checked": 0, "differences": []}
    if declared in CANNOT_JUDGE:
        return {"outcome": "inconclusive", "detail": "W5 could not judge the call: %s%s" % (
            declared, ": %s" % o["note"] if o.get("note") else ""), "checked": 0, "differences": []}
    if declared in NO_CLAIM:
        return {"outcome": "not_checked", "detail": "W5 makes no claim: %s%s" % (
            declared, ": %s" % o["note"] if o.get("note") else ""), "checked": 0, "differences": []}
    if not declared and not n and not o.get("absent"):
        return {"outcome": "not_checked", "detail": "W5 checks nothing of this call (consistency-only)%s" % (
            ": %s" % o["note"] if o.get("note") else ""), "checked": 0, "differences": []}
    code = error_code(other)
    if code in FIT_ERRORS:
        return {"outcome": "not_comparable", "detail": "%s answered %s: a fitting error, about which W5 makes no "
                                                       "claim (FORMAT v2 F3)" % ("the way", code),
                "checked": 0, "differences": []}
    if code is not None and declared != "expected_error":
        return {"outcome": "disagree", "detail": "the error %s where W5 answers (%d checked pointers)" % (code, n),
                "checked": n, "differences": [{"pointer": "/error", "clause": "W5 derives an answer", "oracle":
                                               "an answer", "other": code}]}
    diffs = differences(record_w5, other)
    if declared == "expected_error":
        detail = "an expected error%s" % (" (%s)" % o["expected_error"] if o.get("expected_error") else "")
    else:
        detail = "%d checked pointer(s)" % n
    if diffs:
        d = diffs[0]
        detail = "%s (%s)%s" % (d["pointer"], d.get("clause"), " and %d more" % (len(diffs) - 1) if diffs[1:] else "")
    return {"outcome": "disagree" if diffs else "agree", "detail": detail, "checked": n, "differences": diffs}


# ---------------------------------------------------------------------------------------------- shared derivations
def twin_difference(full_answer: dict, mgt_text: str, seed_index: int = 0) -> Optional[str]:
    """The DESIGN-v8 §2 twin pre-check: the graphlet body against the full twin: the seed (S validated_seed_id and
    the number of seed labels), the label dictionary (L records: name, kind, column, seq_id), and per arm
    (checks.graphlet_arms) status, segments (id, parents, from_bp, length_bp) and paths (id, segments, length_bp,
    end_reasons, path_reason). None when they describe one walk; else the first difference."""
    res = full_answer["results"][seed_index] if isinstance(full_answer.get("results"), list) else full_answer
    try:
        arms = checks.graphlet_arms(mgt_text)
    except (checks.Missing, ValueError, KeyError, IndexError) as e:   # checks.graphlet_arms' refusals
        return "the MGT text does not read: %s" % e
    # the seed and the permitted labels: a time-stopped derivation (SPEC §7.0 D3) can read another number of k-mers in
    # the other execution and walk another label set with the same trie (VO11)
    S = next((x.split(" ", 5) for x in mgt_text.split("\n") if x.startswith("S ")), None)
    seed = res.get("seed") or {}
    if S and len(S) > 4 and "validated_seed_id" in seed and S[1] != seed["validated_seed_id"]:
        return "seed: validated_seed_id full %r, graphlet %r" % (seed["validated_seed_id"], S[1])
    if S and len(S) > 4 and S[4].isdigit() and isinstance(seed.get("labels"), list) and \
            res.get("label_mode") == "constrain" and int(S[4]) != len(seed["labels"]):
        return "seed: %d seed labels in full, %s in the graphlet's S record" % (len(seed["labels"]), S[4])
    try:
        ld = checks.graphlet_label_dict(mgt_text)
    except ValueError as e:
        return "the L records do not read: %s" % e
    fld = res.get("label_dict") or []
    keys = ("name", "kind", "column", "seq_id")
    if len(ld) != len(fld) or any(any(k in fe and fe[k] != ge.get(k) for k in keys) for ge, fe in zip(ld, fld)):
        return "label_dict: %d labels in full, %d L records, or another label" % (len(fld), len(ld))
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


def _drop_timing(o: Any) -> Any:
    if isinstance(o, dict):
        return {k: _drop_timing(v) for k, v in o.items() if k != "timing"}
    if isinstance(o, list):
        return [_drop_timing(x) for x in o]
    return o


def normalise_full_result(result: dict) -> dict:
    """A `detail: full` result as two executions of one walk state it alike (the suite's normalisation, not the
    library's): `timing` removed everywhere and segment events sorted within equal at_bp, needed_budgets sorted
    (DESIGN §9 T37; SPEC §7.1 "events at one at_bp are unordered", §7.2 needed_budgets a histogram); growth
    label_ends null read as {} (DESIGN §2.5); annotation's direct_reads and keys_mapped dropped (SPEC §8.3: they
    count lookahead work, which differs between executions)."""
    r = _drop_timing(copy.deepcopy(result))
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


_TIME_KNOBS = ("bounds.time_budget_ms", "attempt_id")   # attempt_id: an attempt's elapsed ms (SPEC §7.0)
_TIME_RESOURCES = ("time", "attempt_deadline", "cancelled")


def _time_measured(result: dict) -> List[str]:
    """Pointers (relative to a full result) of the values a clock measured: a time or attempt limitation's observed,
    a time resource stop's used, remaining and message (its requested and effective budget are settings, kept: SPEC
    §7.0). Two executions of one walk measure different times (SPEC §6.8: a stop by time is not deterministic;
    §7.2 timing), so W5 checks none of them against its twin."""
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
        out += ["/resource_stop/%s" % k for k in ("used", "remaining", "message") if k in rs]
    return out


_MEMORY_KNOB = "bounds.max_memory_mb"


def _memory_gap(d: _Doc, ctx: dict) -> Optional[str]:
    """The twin ran without the graphlet request's memory budget (VERIFY-D VD4: v8 sends the full twin without
    bounds.max_memory_mb, SPEC §7.1 R21 (4)): the difference, for oracle.inputs, else None. The request's side is
    the graphlet envelope's normalised strategy (context strategy: an omitted budget a server maximum set is in it,
    SPEC §10.3 clamped); the twin's is its envelope's strategy (when it has one) and its result's limitations, which
    must state no memory_bound_soft either (SPEC §7.0: stated exactly when the knob is set). The graphlet answer then
    states the seed-level memory_bound_soft the twin does not: W5 reads its seed-level entries from the graphlet's
    own K records (_seed_limitations)."""
    budget = ((ctx.get("strategy") or {}).get("bounds") or {}).get("max_memory_mb")
    if budget is None:
        return None
    st = d.envelope.get("strategy")
    if isinstance(st, dict) and (st.get("bounds") or {}).get("max_memory_mb") is not None:
        return None
    if any(isinstance(e, dict) and (e.get("kind") == "memory_bound_soft" or e.get("knob") == _MEMORY_KNOB)
           for e in d.result.get("limitations") or []):
        return None
    how = ("its strategy has no bounds.max_memory_mb and " if isinstance(st, dict) else "") + \
        "its result states no memory_bound_soft"
    return "twin without the memory budget: the graphlet request has bounds.max_memory_mb %s, the full twin none " \
           "(%s); the seed-level limitations are the graphlet's own K records (SPEC §7.5.2 K; %s)" % (
               budget, how, C_MEMORY_TWIN)


_K_INT = re.compile(r"(0|[1-9][0-9]*)\Z")
_K_PCT = re.compile(r"%([0-9A-F]{2})")
_K_TEXT = {"%25": "%", "%0A": "\n", "%0D": "\r"}


def _k_value(tok: str) -> Any:
    """A typed K value (SPEC §7.5.1 typed values; §7.5.4 K: "unlimited" is u, every other string s:, an integer i:,
    a real f:): i:<integer>, f:<float> (positional digits, inf), s:<string> (%XX for the bytes outside 0x21-0x7E and
    for % and ,; UTF-8), u."""
    if tok == "u":
        return "unlimited"
    kind, _, v = tok.partition(":")
    if kind == "i" and _K_INT.match(v):
        return int(v)
    if kind == "f" and v:
        return float(v)
    if kind == "s":
        raw = bytearray()
        i = 0
        while i < len(v):
            m = _K_PCT.match(v, i)
            if m:
                raw.append(int(m.group(1), 16))
                i += 3
            else:
                raw += v[i].encode("ascii")
                i += 1
        return raw.decode("utf-8")
    raise ValueError("K value %r: i:, f:, s: or u (SPEC §7.5.1)" % tok)


def k_entries(mgt_text: str) -> List[Tuple[str, dict]]:
    """The graphlet's K records as (arm: "*" seed level, "l", "r"; the JSON entry), in the body's order (SPEC §7.5.2
    K: "one entry of limitations (§7.0), the JSON entry itself", seed-level entries first, then the left arm's, then
    the right arm's, each list in its JSON order): kind, knob, limit, observed, complete_to_bp (`*` absent), the extra
    fields, effect (free text: %25, %0A, %0D, SPEC §7.5.1)."""
    out = []
    for line in mgt_text.split("\n"):
        if not line.startswith("K "):
            continue
        f = line.split(" ", 8)
        if len(f) != 9:
            raise ValueError("K record with %d fields: %r" % (len(f), line[:80]))
        e: dict = {"kind": f[2], "knob": f[3]}
        for key, tok in (("limit", f[4]), ("observed", f[5])):
            if tok != "*":
                e[key] = _k_value(tok)
        if f[6] != "*":
            if not _K_INT.match(f[6]):
                raise ValueError("K complete_to_bp %r" % f[6])
            e["complete_to_bp"] = int(f[6])
        if f[7] != ".":
            for item in f[7].split(","):
                name, eq, tok = item.partition("=")
                if not eq or not name:
                    raise ValueError("K extra %r: name=VALUE" % item)
                e[name] = _k_value(tok)
        e["effect"] = re.sub(r"%25|%0A|%0D", lambda m: _K_TEXT[m.group(0)], f[8])
        out.append((f[1], e))
    return out


def _seed_limitations(d: _Doc, ctx: dict) -> Optional[List[dict]]:
    """Under a twin without the memory budget (_memory_gap), the seed-level limitations as the graphlet states them:
    its K records with arm * (SPEC §7.5.2 K), memory_bound_soft included; else None (the twin's are the graphlet's)."""
    if not _memory_gap(d, ctx):
        return None
    return [e for arm, e in k_entries(d.mgt) if arm == "*"]


C_SEED_K = "SPEC §7.5.2 K (arm *): the graphlet's own seed-level entries, each the JSON entry itself, in their " \
           "order (the twin ran without the memory budget: " + C_MEMORY_TWIN + ")"


def _evidence(d: _Doc, sides: Sequence[str]) -> Tuple[dict, Dict[str, str], List[str]]:
    """The evidence block of every local answer (DESIGN §6: complete_to_bp, exact, support, reconverge, scope,
    outcome, limitations), for the arms the answer covers."""
    ev: dict = {"complete_to_bp": {s: d.arms[s]["complete_to_bp"] for s in sides},
                "scope": {s: d.arms[s]["completeness_scope"] for s in sides},
                "outcome": d.result["outcome"]}
    if d.h(7) in _SUPPORT:
        ev["support"] = _SUPPORT[d.h(7)]
    if d.h(8) in _RECONVERGE:
        ev["reconverge"] = _RECONVERGE[d.h(8)]
    clauses = {"/evidence/complete_to_bp": C_CTB, "/evidence/scope": C_SCOPE, "/evidence/outcome": C_OUTCOME,
               "/evidence/support": "SPEC §7.5.2 H support k|t (kmer, trace)",
               "/evidence/reconverge": "SPEC §7.5.2 H reconverge m|k (merge, keep)"}
    cons = ["/evidence/coordinates", "/evidence/view", "/evidence/qualified", "/evidence/limitations"]
    # DESIGN §6 states a necessary condition only ("exact requires complete label evidence and no cut lists"): where
    # it fails, exact is false; where it holds, exact is the tool's (consistency-only)
    cut = any((d.arms[s].get("labels_per_node") or {}).get("nodes_truncated", 0) > 0 for s in sides)
    if (d.result.get("outcome") or {}).get("label_evidence") != "complete" or cut:
        ev["exact"] = False
        clauses["/evidence/exact"] = "DESIGN §6: exact requires complete label evidence and no cut lists " \
                                     "(outcome.label_evidence, labels_per_node.nodes_truncated)"
    else:
        cons.append("/evidence/exact")
    # evidence.limitations is the tool's digest of the stated limitations (finding F1: it leaves out a scope entry
    # that observed 0 merges and adds the resource stop): DESIGN §6 names the key without defining it
    return ev, clauses, cons


def _cut_rows(ans: dict, ctx: dict, consistency: List[str], key: str = "rows") -> List[str]:
    """Whole lists unless the op read fewer pages: then the reference way's row count, or no rows at all. -> the
    partial list's pointer (for _record), or []."""
    if ctx.get("paging", "all") == "all":
        return []
    if isinstance(ctx.get("rows"), int):
        ans[key] = ans[key][:ctx["rows"]]
    else:
        ans[key] = []
        consistency.append("/%s/*" % key)
    return ["/" + key]


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
         "index": {"ns": d.h(13), "verifiable": d.h(14) is not None}, "arms": {}}
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
    arm's, each in its JSON order (SPEC §7.5.2 K), and the resource stop (SPEC §7.0).
    A twin run without the graphlet request's memory budget (_memory_gap) states no memory_bound_soft, which the
    answer holds (SPEC §7.0): the seed-level entries are then the graphlet's own K records (SPEC §7.5.2 K, the JSON
    entries in their order, memory_bound_soft included: _seed_limitations), checked as they are."""
    seed_k = _seed_limitations(d, ctx)
    seed_level = list(d.result.get("limitations") or []) if seed_k is None else seed_k
    # (the entry's pointer in the twin's result or None, the row): seed-level entries carry no arm
    src = [("/limitations/%d" % j if seed_k is None else None, dict(e)) for j, e in enumerate(seed_level)] \
        if side is None else []
    for s in d.sides:
        if side is None or s == side:
            src += [("/arms/%s/limitations/%d" % (s, j), dict(e, arm=s))
                    for j, e in enumerate(d.arms[s].get("limitations") or [])]
    rows = [r for _, r in src]
    timed = set(_time_measured(d.result))
    cons_t = []
    for i, (base, row) in enumerate(src):
        if base is not None and base + "/observed" in timed:
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
    if seed_k is not None and side is None:
        clauses.update({"/rows/%d/*" % i: C_SEED_K for i in range(len(seed_k))})
    clauses.update({"/total": "SPEC §7.0: every stated limitation (seed level and the arms asked)" +
                              ("; the seed level from the graphlet's K records (%s)" % C_SEED_K
                               if seed_k is not None else ""),
                    "/rows/*/arm": "SPEC §7.5.2 K arm (an arm's entry names its arm)",
                    "/rows/*/*": "SPEC §7.0 limitations: the stated entry (§7.5.2 K: the JSON entry itself)"})
    cons = ["/rows/*/informational", "/next_cursor"] + cons_t
    if side is not None and seed_level:
        # whether a seed-level entry belongs to an arm's listing is the tool's choice
        ans.pop("total")
        ans["rows"] = []
        cons += ["/total", "/rows/*"]
    partial = _cut_rows(ans, ctx, cons)
    return _record(ans, clauses, cons, partial=partial)


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
             "length": "DESIGN §5 names by='length' without a key: longest first, ties by path_id (INVENTORY "
                       "§2.6, presentation)",
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
    partial = _cut_rows(ans, ctx, cons) if "rows" in ans else []
    ev, ecl, econ = _evidence(d, [side])
    ans["evidence"] = ev
    clauses.update(ecl)
    # computed with DESIGN §5.1's displayed-evidence formula from engine values (end_labels[].route_bp, runs'
    # from_bp), rather than read from the engine: the report marks them "re-derived from DESIGN" (DESIGN-v8 §2)
    red = ["/rows/*/n_full", "/rows/*/labels_full", "/rows/*/labels_full_truncated", "/filtered/merge_entered"]
    if rank == "support":
        red.append("/rows/*/walk")
    if lab is not None and rc:
        red += ["/total", "/rows/*/walk"]
    return _record(ans, clauses, cons + econ, partial=partial, rederived=red)


def _chain_rows(d: _Doc, side: str, p: dict, cap: Optional[int]) -> List[dict]:
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


def _ends_order(pattern: str) -> dict:
    """How a walk row's `ends` compare: as a set while ENDS_UNORDERED (G3), else sorted by position only (SPEC
    §7.1)."""
    return {"unordered": [pattern]} if ENDS_UNORDERED else {"position_ordered": [pattern]}


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
    rows = _chain_rows(d, side, p, LIST_CAP)
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
    partial = _cut_rows(ans, ctx, cons)
    kept = len(ans["rows"])   # the row-level absent checks of rows W5 still answers stay (VO15)
    absent = {k: v for k, v in absent.items() if not k.startswith("/rows/") or int(k.split("/")[2]) < kept}
    ev, ecl, econ = _evidence(d, [side])
    ans["evidence"] = ev
    clauses.update(ecl)
    return _record(ans, clauses, cons + econ, absent=absent, partial=partial, **_ends_order("/rows/*/ends"))


def _holds_base(a: dict, label: int, at_bp: int, mode: str) -> Optional[bool]:
    """Whether a label is present at outward base |at_bp| of the arm: constrain, one of its runs holds the base
    (SPEC §7.1 runs[], half-open [from_bp, to_bp)); annotate, a recorded label set holds it (SPEC §7.4 label_sets,
    [from_bp, to_bp) covering the nodes entered by steps from_bp + 1 .. to_bp). None where a cut set covers the base
    without the label (a lower bound: SPEC §7.4)."""
    if mode == "constrain":
        return any(r["label"] == label and r["from_bp"] <= at_bp < r["to_bp"] for r in a.get("runs") or [])
    unknown = False
    for sg in a.get("segments") or []:
        for ls in sg.get("label_sets") or []:
            if ls["from_bp"] <= at_bp < ls["to_bp"]:
                if label in (ls.get("labels") or []):
                    return True
                unknown = unknown or bool(ls.get("truncated"))
    return None if unknown else False


def _graphlet_labels(d: _Doc, args: dict, ctx: dict) -> dict:
    if args.get("name") is not None:
        return _label_walks(d, args, ctx)
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
    undecided = False
    for e in d.result.get("label_summary") or []:
        for s in sides:
            v = e.get(s)
            if not isinstance(v, dict) or v["direct_bp"] < min_direct:
                continue
            if at_bp is not None:
                held = _holds_base(d.arms[s], e["label"], at_bp, d.mode)
                if held is None:
                    undecided = True
                if not held:
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
    if undecided:
        # a cut label set covers the base: which labels are present there is a lower bound (SPEC §7.4)
        ans = {}
        cons += ["/total", "/rows/*"]
    else:
        if at_bp is not None:
            clauses["/total"] += "; at_bp: the label is present at the base (a run holds it, SPEC §7.1; annotate: " \
                                 "a recorded label set, SPEC §7.4)"
    partial = _cut_rows(ans, ctx, cons) if "rows" in ans else []
    ev, ecl, econ = _evidence(d, sides)
    ans["evidence"] = ev
    clauses.update(ecl)
    return _record(ans, clauses, cons + econ, partial=partial)


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
    partial = _cut_rows(ans, ctx, cons)
    ev, ecl, econ = _evidence(d, [side])
    ans["evidence"] = ev
    clauses.update(ecl)
    return _record(ans, clauses, cons + econ, partial=partial)


def _graphlet_sequence(d: _Doc, args: dict, ctx: dict) -> dict:
    side = d.arm(args.get("arm"), required=True)
    # typed as the tool types them (library-defined, INVENTORY §1.1, §2.12): integers never booleans
    walk, segment = _int_arg(args, "walk", None, 0), _int_arg(args, "segment", None, 0)
    orientation = args.get("orientation", "natural")
    with_seed = _bool_arg(args, "with_seed", False)
    start, end = _int_arg(args, "start", 0, 0), _int_arg(args, "end", None, 0)
    a = d.arms[side]
    if (walk is None) == (segment is None) or orientation not in ("natural", "walk") or \
            (with_seed and segment is not None) or (end is not None and end < start):
        raise _Refused("walk %r / segment %r / orientation %r / with_seed %r / [%r, %r)" % (
            walk, segment, orientation, with_seed, start, end))
    if walk is not None:
        if not walk < len(a["paths"]):
            raise _Refused("walk %r" % (walk,))
        flanks = d.flanks(side)
        if flanks is None:
            raise _Refused("no bases (sequences: false)", "no_bases", C_NO_BASES)
        natural = d.molecule(side, flanks[walk]) if with_seed else flanks[walk]
    else:
        if not segment < len(a["segments"]):
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
    ev, ecl, econ = _evidence(d, [side])
    page = len(json.dumps(dict(ans, evidence=ev, to=stop), separators=(",", ":"))) + _SEQ_PAGE_EXTRA
    fits = stop - start + page <= ctx.get("sequence_max_bytes", 16384)
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
    """The overlap of a GFA CIGAR `nM`, else None."""
    return int(c[:-1]) if c and c[-1] == "M" and c[:-1].isdigit() else None


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
    refused = _path_refused(args, ctx)
    if refused is not None:
        return refused
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
                           "library-defined: what=walk needs arm and walks, and the code (INVENTORY §2.13); DESIGN §6 "
                           "only says what=walk writes the complete chains")
        side = d.arm(args.get("arm"), required=True)
        ids = _walk_ids(d, side, args.get("walks"))
        walks, absent = [], {}
        for j, w in enumerate(ids):
            p = d.arms[side]["paths"][w]
            rows = _chain_rows(d, side, p, None)
            absent.update(_walk_row_extras(rows, "/walks/%d/segments" % j, None)[1])
            walks.append({"walk": w, "length_bp": p["length_bp"], "segments": rows})
        clauses = {"/arm": C_ARM_ECHO, "/walks/*/walk": C_PATHS, "/walks/*/length_bp": C_LENGTH}
        clauses.update({"/walks/*/segments/*/" + k: v.replace("first 8", "uncut").replace(
            "checked when at most 8 (INVENTORY §2.7)", "uncut") for k, v in _WALK_ROW_CLAUSES.items()})
        return _record({"arm": side, "walks": walks}, clauses,
                       cons + ["/records", "/walks/*/sequence", "/walks/*/segments/*/exact", "/handle", "/evidence"],
                       absent=absent, **_ends_order("/walks/*/segments/*/ends"))
    if fmt == "json":
        res = normalise_full_result(d.result)
        c_json = "SPEC §7.5.6: the full result is reproduced exactly from the text and its summary (the twin's " \
                 "result, normalised: DESIGN §9 T37, §2.5, SPEC §8.3)"
        timed = ["/json" + t for t in _time_measured(res)]
        seed_k = _seed_limitations(d, ctx)
        if seed_k is not None:
            # the twin states no memory_bound_soft (it ran without the budget): the seed-level list is the graphlet's
            # own K records (C_SEED_K), its clock values included (they are the graphlet's statement)
            res["limitations"] = seed_k
            timed = [t for t in timed if not t.startswith("/json/limitations/")]
        for t in timed:   # a clock's value: left out of W5's answer, consistency-only (SPEC §6.8)
            parent, key = t.rsplit("/", 1)
            resolve_pointer({"json": res}, parent).pop(key.replace("~1", "/").replace("~0", "~"), None)
        clauses = {}
        for k in res:
            if k == "limitations" and seed_k is not None:
                clauses["/json/limitations"] = C_SEED_K
            elif k == "arms":
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
                # a GFA 1.0 file can carry `*` sequences: whether the tool writes one or refuses is its own choice
                return _record({}, {}, cons + ["/records", "/gfa"], note="gfa of a retrieval without bases: "
                                                                         "consistency-only (GFA allows `*`)")
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
            raise _Refused("fasta of a retrieval without bases", None, C_NO_BASES_EXPORT)
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
    refused = _path_refused(args, ctx)
    if refused is not None:
        return refused
    return _record({"body_digest": body_digest(d.mgt), "arms_digest": arms_digest(d.mgt)},
                   {"/body_digest": "SPEC §7.5.2: the saved file is H, J, the body unchanged (FORMAT-v8 §3 body_digest)",
                    "/arms_digest": "the saved file reads as the original (checks.graphlet_arms, SPEC §7.5)"},
                   ["/path", "/bytes", "/parsed", "/validation", "/fields_cut"], inputs=("mgt",))


def _graphlet_load(d: _Doc, args: dict, ctx: dict) -> dict:
    refused = _path_refused(args, ctx, spool=True)
    if refused is not None:
        return refused
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


# ---------------------------------------------------------------------------------------------- phase B
# SPEC §7.5.1: a label end's text qualifier, the second letter of its code (Lw, a silent end, and m, a merge's
# close, are end tokens of R, not qualifiers)
_QUALIFIER = {"Rm": "minority", "Rb": "below_min_labels", "Rs": "split_limit", "Dh": "hairpin", "Ls": "superseded",
              "Lx": "switch_sources"}
C_EVIDENCE = "DESIGN §5.1 displayed-path evidence (pseudocode, re-derived): the latest merge on the first-parent " \
             "chain to the run's anchor whose first partition (labels_via_parent[0], SPEC §7.1) lacks the lineage's " \
             "incoming label (R prev_run), at least the run's from_bp"
C_RUN = "SPEC §7.5.2 R / §7.1 runs[] (the engine's run, run id order)"


def _constrain_runs(d: _Doc, side: str) -> bool:
    """Whether the arm's runs can be read lineage by lineage: constrain mode, R records aligned with runs[]."""
    return d.mode == "constrain" and d.aligned(side)


def _selected_labels(d: _Doc, sel: Any) -> Optional[set]:
    """DESIGN §5 Graphlet.claims(labels=): one selector or a list of them; None = every label."""
    if sel is None:
        return None
    sels = sel if isinstance(sel, list) else [sel]
    return {d.resolve(x) for x in sels}


def _graphlet_claims(d: _Doc, args: dict, ctx: dict) -> dict:
    """The arm's claims (DESIGN §5 Graphlet.claims; SPEC §6.9: one claim per RUN END): per run of the twin, in run id
    order, its label, interval, route stamp, entry, switch source and cost, loss and branches (the engine's run),
    the end reason (the run's end_reason; the qualifier from R's end code), the path a leaf run ends on, and its
    displayed evidence (DESIGN §5.1, re-derived). Filters: labels (DESIGN §5), route_consistent (DESIGN §6: a claim
    whose SPEC §7.1 merge stamp is above 0, T at a leaf and R inside, is left out and counted merge_entered: the
    reading phase A's walks filter uses, F2, VO-B2; where R's and the displayed stamp disagree the selection is
    consistency-only), min_bp (its displayed support to_bp - evidence_from at least: INVENTORY §2.11). Annotate claims (DESIGN §5.2.1's witness rule) are consistency-only."""
    side = d.arm(args.get("arm"), required=True)
    min_bp = _int_arg(args, "min_bp", 0, 0)
    rc = _bool_arg(args, "route_consistent", True)
    if args.get("labels") == []:
        # DESIGN §5 defines labels=None as every label; what an empty selection means is the tool's
        return _record({}, {}, ["/*"], note="labels=[]: library-defined, consistency-only")
    labs = _selected_labels(d, args.get("labels"))
    cons = ["/rows/*/end_class", "/rows/*/support", "/rows/*/exact", "/rows/*/coordinates",
            "/rows/*/coordinates_total", "/rows/*/coordinates_truncated", "/rows/*/coordinates_kind",
            "/rows/*/coordinates_lower_bound", "/hint", "/caveat", "/next_cursor", "/coordinates"]
    ev, ecl, econ = _evidence(d, [side])
    if not _constrain_runs(d, side):
        return _record({"arm": side, "evidence": ev}, dict(ecl, **{"/arm": C_ARM_ECHO}),
                       cons + econ + ["/total", "/rows", "/filtered", "/rows/*/kind"],
                       note="annotate claims follow DESIGN §5.2.1's witness rule: consistency-only")
    a = d.arms[side]
    leaf_len = {p["segments"][-1]: p for p in a["paths"]}
    # a run's R stamp is the earliest non-first-parent merge on its lineage's route, a leaf label's T stamp the latest
    # (SPEC §7.1); DESIGN does not say which a claim carries (finding PB2): W5 checks route_bp where R's stamp equals
    # the displayed one (T at a leaf, the evidence's merge part inside), and leaves it open where they differ
    leaf_t = {e["run"]: e.get("route_bp", 0) for p in a["paths"] for e in p.get("end_labels") or []
              if isinstance(e.get("run"), int)}
    rows, filtered, absent = [], {}, {}
    open_stamp = []
    undecided_at: Optional[int] = None
    for i, run in enumerate(d.runs(side)):
        if labs is not None and run["label"] not in labs:
            continue
        ef = d.evidence_from(side, i)
        if ef is None:
            return _record({"arm": side, "evidence": ev}, dict(ecl, **{"/arm": C_ARM_ECHO}),
                           cons + econ + ["/total", "/rows", "/filtered", "/rows/*/kind"],
                           note="a merge on a claim's chain states no partitions: consistency-only")
        if rc:
            # one route filter for walks and claims (VO-B2, phase A F2): merge-entered where the claim's SPEC §7.1
            # stamp is above 0 (a leaf label's T stamp, else the run's R stamp); where R's stamp and the displayed
            # one (the evidence's merge part) disagree on that, no clause decides: the selection is open
            stamp = leaf_t[i] if i in leaf_t else run.get("route_bp", 0)
            other = leaf_t[i] if i in leaf_t else d.evidence_from(side, i, merge_part=True)
            if (stamp > 0) != ((other or 0) > 0) or (i in leaf_t and (run.get("route_bp", 0) > 0) != (stamp > 0)):
                undecided_at = len(rows)   # this claim and every later row: open (R2B-3)
                break
            if stamp > 0:
                filtered["merge_entered"] = filtered.get("merge_entered", 0) + 1
                continue
        if run["to_bp"] - ef < min_bp:
            continue
        code = d.end_code(side, i)
        boundary = run["from_bp"] == run["to_bp"]
        row: dict = {"label": d.label(run["label"]), "from_bp": run["from_bp"], "to_bp": run["to_bp"],
                     # a boundary claim has no bases, so no displayed interval: null (DESIGN §5.1; its clipped-merge
                     # vector writes null wherever nothing is displayed; finding PB1)
                     "route_bp": run.get("route_bp", 0), "evidence_from": None if boundary else ef,
                     "entered_by": run["entered_by"],
                     "from_label": d.label(run["from"]) if run.get("from") is not None else None,
                     "cost": run.get("cost"), "loss": run.get("loss"), "branches": run.get("branches"),
                     "qualifier": _QUALIFIER.get(code)}
        if "end_reason" in run:
            row["reason"] = run["end_reason"]
        p = leaf_len.get(run["segment"])
        if p is not None and p["length_bp"] == run["to_bp"] and code != "m":
            row["path_id"] = p["id"]
        if boundary:
            row["kind"] = "boundary"
        shown = leaf_t[i] if i in leaf_t else d.evidence_from(side, i, merge_part=True)
        if shown != row["route_bp"]:
            del row["route_bp"]
            open_stamp.append(len(rows))
        rows.append(row)
    cons += ["/rows/%d/route_bp" % j for j in open_stamp]
    for j, row in enumerate(rows):
        if "path_id" not in row:
            cons.append("/rows/%d/path_id" % j)
        if "reason" not in row:
            cons.append("/rows/%d/reason" % j)
        if "kind" not in row:   # DESIGN names the other kinds without defining them
            cons.append("/rows/%d/kind" % j)
    ans: dict = {"arm": side, "total": len(rows), "rows": rows, "evidence": ev}
    note = None
    if undecided_at is not None:
        # a claim whose R stamp and displayed stamp disagree on a non-first-parent entry (PB2, VO-B2): whether the
        # filter keeps it is open, so it and every later row, the total and the counts are consistency-only; the
        # rows before it stay checked (R2B-3)
        del ans["total"]
        cons += ["/total", "/filtered"]
        note = "route_consistent: rows from %d on are consistency-only (an undecided claim)" % undecided_at
    elif filtered and not min_bp:
        ans["filtered"] = filtered
    elif min_bp:
        # whether the route filter counts before or after min_bp is the tool's order (finding PB3)
        cons.append("/filtered")
    clauses = {"/arm": C_ARM_ECHO,
               "/total": "SPEC §6.9: one claim per run end, the filters kept (DESIGN §5 labels; route_consistent and "
                         "min_bp on the displayed evidence, DESIGN §5.1, §6)",
               "/filtered/*": "DESIGN §6: a filter never hides rows silently (merge_entered: the claim's SPEC §7.1 "
                              "stamp is above 0, a leaf label's T stamp or the run's R stamp: phase A F2's reading)",
               "/rows/*/label": "the run's label (" + C_RUN + "); " + C_LABEL,
               "/rows/*/from_bp": C_RUN, "/rows/*/to_bp": C_RUN,
               "/rows/*/route_bp": "SPEC §7.1, §7.5.2 R route_bp: the run's earliest non-first-parent merge stamp "
                                   "(where the displayed stamp differs, the T stamp at a leaf or the evidence's merge "
                                   "part inside, which one a claim carries is open: consistency-only, finding PB2)",
               "/rows/*/evidence_from": C_EVIDENCE + "; null for a boundary claim (no bases, DESIGN §5.1)",
               "/rows/*/entered_by": "SPEC §7.1 runs[].entered_by (seed | switch)",
               "/rows/*/from_label": "SPEC §7.1 runs[].from (the switch source; null when entered by seed); " + C_LABEL,
               "/rows/*/cost": "SPEC §7.5.2 R from_label:cost (null when entered by seed)",
               "/rows/*/loss": "SPEC §7.1 runs[].loss (the lineage's at the end, §4)",
               "/rows/*/branches": "SPEC §7.1 runs[].branches (the lineage's at the end)",
               "/rows/*/reason": "SPEC §7.1 runs[].end_reason (the enum; a merge's close states none)",
               "/rows/*/qualifier": "SPEC §7.5.1: the end code's text qualifier (Rm, Rb, Rs, Dh, Ls, Lx), else null",
               "/rows/*/path_id": "SPEC §7.5.6: a run anchored at a leaf with to_bp = the path length is that path's "
                                  "end label (the path a claim of an inner run is shown on: no clause)",
               "/rows/*/kind": "DESIGN §5.1: a run of zero length is a boundary claim"}
    clauses.update(ecl)
    red = ["/rows/*/evidence_from"]
    if min_bp:   # min_bp selects on the displayed evidence (DESIGN §5.1, re-derived); the route filter on stamps
        red += ["/total", "/rows/*/label", "/rows/*/from_bp", "/rows/*/to_bp"]
    partial = _cut_rows(ans, ctx, cons)
    if undecided_at is not None:
        partial = ["/rows"]
    return _record(ans, clauses, cons + econ, partial=partial, rederived=red, note=note)


def _merged_into(d: _Doc, side: str, run: dict) -> Optional[int]:
    """The merge that closed a run (R end m): the segment created at to_bp with the run's anchor among its parents
    (SPEC §7.5.2 R m; DESIGN §4: every m run's anchor is a parent of the segment created at to_bp)."""
    hits = [sg["id"] for sg in d.arms[side]["segments"]
            if len(sg.get("parents") or []) > 1 and sg["from_bp"] == run["to_bp"] and run["segment"] in sg["parents"]]
    return hits[0] if len(hits) == 1 else None


def _label_walks(d: _Doc, args: dict, ctx: dict) -> dict:
    """graphlet_labels(name=): the named label's runs (DESIGN §5 label_walks, §6 "that label's runs/routes incl.
    merge routes"), per arm (both, left first, when arm is None), in run id order: from_bp, to_bp and end (the
    engine's run, SPEC §7.1 runs[]; the end named as SPEC §7.5.6 names a label end), evidence_from (DESIGN §5.1,
    re-derived, also for a zero-length run: PB4), merged_into (the merge that closed it, SPEC §7.5.2 R m) and route
    (DESIGN §5 routes, re-derived); walks_below is consistency-only (DESIGN names leaves_below without defining
    it: VO-B3). rank, min_direct_bp and at_bp do not apply (INVENTORY §2.9). Annotate: the
    witness route of DESIGN §5.2.1, consistency-only."""
    lab = d.resolve(args["name"])
    side = d.arm(args.get("arm"), required=False)
    sides = [side] if side else d.sides
    if _bool_arg(args, "coordinates", False) and not isinstance(d.result.get("coordinates"), dict):
        raise _Refused("coordinates=true on a retrieval without record coordinates", "no_coordinates",
                       "SPEC §7.1 record coordinates (null, coordinates_reason); the code no_coordinates is "
                       "library-defined (INVENTORY §2.9)")
    ev, ecl, econ = _evidence(d, sides)
    base = {"name": d.labels[lab]["name"], "ref": d.labels[lab]["ref"]}
    c_base = {"/name": "the label named (" + C_LABEL + ")", "/ref": "the label named (" + C_LABEL + ")"}
    # walks_below: DESIGN §5 gives only the name leaves_below; leaves reachable from the anchor and leaves whose
    # displayed chain passes it differ (on 184 of 3,730 fixture runs), and no clause chooses (VO-B3, a DESIGN gap)
    cons = ["/rows/*/coordinates", "/rows/*/coordinates_total", "/rows/*/coordinates_truncated",
            "/rows/*/coordinates_kind", "/rows/*/coordinates_lower_bound", "/next_cursor", "/coordinates",
            "/rows/*/walks_below"]
    if not all(_constrain_runs(d, s) for s in sides):
        return _record(dict(base, evidence=ev), dict(ecl, **c_base), cons + econ + ["/total", "/rows"],
                       note="annotate: the witness route of DESIGN §5.2.1, consistency-only")
    rows, red_rows = [], []
    for s in sides:
        for i, run in enumerate(d.runs(s)):
            if run["label"] != lab:
                continue
            ef = d.evidence_from(s, i)
            route = d.route(s, i)
            code = d.end_code(s, i)
            # a run, not a claim: DESIGN §5.1's formula, also for a zero-length run (the boundary-claim rule is the
            # claims'; finding PB4)
            row: dict = {"arm": s, "from_bp": run["from_bp"], "to_bp": run["to_bp"], "evidence_from": ef,
                         "merged_into": _merged_into(d, s, run) if code == "m" else None}
            if "end_reason" in run and code not in ("Lw", "m"):
                # SPEC §7.5.6: a label end's reason is the qualifier's text or the code's enum
                row["end"] = _QUALIFIER.get(code) or run["end_reason"]
            else:   # a silent switch-source end (Lw) and a merge's close have no label_end event: no SPEC word
                cons.append("/rows/%d/end" % len(rows))
            if route is not None:
                row["route"] = route
            else:
                cons.append("/rows/%d/route" % len(rows))
            if ef is None:
                del row["evidence_from"]
                cons.append("/rows/%d/evidence_from" % len(rows))
            rows.append(row)
    ans: dict = dict(base, total=len(rows), rows=rows, evidence=ev)
    clauses = dict(c_base, **{
        "/total": "the named label's runs on the arms asked (SPEC §7.5.2 R, §7.1 runs[])",
        "/rows/*/arm": "the arm of the run (" + C_ARM_ECHO + "), left before right (SPEC §7.5.2 document order)",
        "/rows/*/from_bp": C_RUN, "/rows/*/to_bp": C_RUN,
        "/rows/*/end": "SPEC §7.5.6: the label end's reason, the qualifier's text (R's end code, SPEC §7.5.1) or "
                       "the enum (runs[].end_reason); Lw and m have no label_end event: c/o",
        "/rows/*/evidence_from": C_EVIDENCE + " (a run, not a claim: the formula also for a zero-length run, PB4)",
        "/rows/*/merged_into": "SPEC §7.5.2 R m and DESIGN §4: the merge segment created at to_bp with the run's "
                               "anchor among its parents; null when no merge closed the run",
        "/rows/*/route": "DESIGN §5 routes (re-derived): root -> the run's anchor, at each merge the parent whose "
                         "partition (labels_via_parent, SPEC §7.1) holds the lineage's incoming label"})
    clauses.update(ecl)
    partial = _cut_rows(ans, ctx, cons)
    return _record(ans, clauses, cons + econ, partial=partial, rederived=["/rows/*/evidence_from", "/rows/*/route"])


def _support_sets(d: _Doc, side: str, p: dict) -> Optional[List[set]]:
    """Per outward base of walk |p| (0 .. length-1), the labels with displayed support there (DESIGN §5 support
    profile, kind displayed; re-derived from DESIGN §5.1): each run's displayed evidence, derived at its anchored
    endpoint with p = the first-parent chain root -> its anchor ([evidence_from, to_bp), "the displayed bases before
    s are not this lineage's"), kept where that chain is the walk (the walk's segment at the base lies on it). The
    library's presence profile (a run's route passing the walk's segment) is library defect G4 (VERIFY-O-B VO-B1).
    None where a merge states no partitions."""
    seg_at: List[int] = []
    for sid in p["segments"]:
        seg_at += [sid] * d.segment(side, sid)["length_bp"]
    sets: List[set] = [set() for _ in range(p["length_bp"])]
    for i, run in enumerate(d.runs(side)):
        if run["from_bp"] == run["to_bp"]:
            continue
        ef = d.evidence_from(side, i)
        if ef is None:
            return None
        chain = set(d.chain(side, run["segment"]))
        for b in range(ef, min(run["to_bp"], p["length_bp"])):
            if seg_at[b] in chain:
                sets[b].add(run["label"])
    return sets


def _recorded_sets(d: _Doc, side: str, p: dict) -> Optional[List[Tuple[Tuple[int, ...], int, bool]]]:
    """annotate: per outward base of walk |p|, the label set recorded at the walk's node there (SPEC §7.4 label_sets:
    [from_bp, to_bp) covers the nodes entered by steps from_bp + 1 .. to_bp), as (the recorded list, its true count,
    cut). None where a base has no recorded set."""
    per: List[Optional[Tuple[Tuple[int, ...], int, bool]]] = [None] * p["length_bp"]
    for sid in p["segments"]:
        for ls in d.segment(side, sid).get("label_sets") or []:
            for b in range(ls["from_bp"], min(ls["to_bp"], p["length_bp"])):
                per[b] = (tuple(ls.get("labels") or []), ls.get("labels_total", len(ls.get("labels") or [])),
                          bool(ls.get("truncated")))
    return None if any(x is None for x in per) else per  # type: ignore[return-value]


def _graphlet_support(d: _Doc, args: dict, ctx: dict) -> dict:
    """graphlet_support: the support profile along a walk (DESIGN §5 support_profile(kind='displayed'), §6 "support
    runs with who left / joined"): maximal stretches [from_bp, to_bp) of one set of labels with displayed support
    (re-derived from DESIGN §5.1, _support_sets; a stretch no label supports is an empty set, as SPEC §6.11 records
    one), each with
    n_labels and, after the first, the labels that joined and left (ascending label id, first 8). Consistency-only:
    the first row's lists, why, exact, n_labels_total on cut lists, every row with step= (its merging rule is the
    tool's). In annotate mode the profile is the labels the engine recorded at the walk's nodes (SPEC §7.4
    label_sets: engine-computed, not re-derived); who joined or left is not read across a cut list."""
    side = d.arm(args.get("arm"), required=True)
    w = args.get("walk")
    paths = d.arms[side]["paths"]
    if not isinstance(w, int) or isinstance(w, bool) or not 0 <= w < len(paths):
        raise _Refused("walk %r of %d" % (w, len(paths)))
    step = _int_arg(args, "step", None, 0)
    ev, ecl, econ = _evidence(d, [side])
    cons = ["/rows/*/why", "/rows/*/exact", "/rows/0/joined", "/rows/0/left", "/next_cursor", "/coordinates"]
    base = {"arm": side, "walk": w}
    c_base = {"/arm": C_ARM_ECHO, "/walk": C_PATHS}
    annotate = d.mode == "annotate"
    if annotate:
        # SPEC §7.4: the labels recorded at every node of the walk (engine-computed), maximal runs of one (capped
        # list, count); a cut list states its true count
        rec_ = _recorded_sets(d, side, paths[w])
        sets = None if rec_ is None else [(set(x[0]), x[1], x[2]) for x in rec_]
    else:
        lin = _support_sets(d, side, paths[w]) if _constrain_runs(d, side) else None
        sets = None if lin is None else [(x, len(x), False) for x in lin]
    if sets is None or step is not None:
        why = "step=: the tool's merging rule" if sets is not None else "no recorded set or partition to read"
        return _record(dict(base, evidence=ev), dict(ecl, **c_base), cons + econ + ["/total", "/rows"],
                       note="%s: consistency-only" % why)
    rows: List[dict] = []
    for b, (cur, total, cut_list) in enumerate(sets):
        if rows and (rows[-1]["_set"], rows[-1]["_total"]) == (cur, total):
            rows[-1]["to_bp"] = b + 1
            rows[-1]["_cut"] = rows[-1]["_cut"] or cut_list
        else:
            rows.append({"from_bp": b, "to_bp": b + 1, "_set": cur, "_total": total, "_cut": cut_list})
    out, prev = [], None
    cut = not annotate and (d.arms[side].get("labels_per_node") or {}).get("nodes_truncated", 0) > 0
    for j, r in enumerate(rows):
        row = {"from_bp": r["from_bp"], "to_bp": r["to_bp"], "n_labels": len(r["_set"])}
        if not cut:
            row["n_labels_total"] = r["_total"]
        if prev is not None:
            if r["_cut"] or prev["_cut"]:   # a cut list: who joined or left cannot be read from it
                cons += ["/rows/%d/joined" % j, "/rows/%d/left" % j]
            else:
                for key, ids in (("joined", sorted(r["_set"] - prev["_set"])), ("left", sorted(prev["_set"] - r["_set"]))):
                    if len(ids) > LIST_CAP:
                        cons.append("/rows/%d/%s" % (j, key))
                    elif ids:
                        row[key] = [d.label(x) for x in ids]
        prev = r
        out.append(row)
    ans: dict = dict(base, total=len(out), rows=out, evidence=ev)
    c_sup = "DESIGN §5 support_profile (displayed; re-derived from DESIGN §5.1): each run's [evidence_from, to_bp) " \
            "derived at its anchor (p = the first-parent chain to it), kept where that chain is the walk; maximal " \
            "stretches of one label set (an empty set too, as SPEC §6.11 records one)"
    if annotate:
        c_sup = "SPEC §7.4 label_sets: the labels recorded at the walk's nodes (annotate), maximal runs of one " \
                "(capped list, count)"
    clauses = dict(c_base, **{"/total": c_sup, "/rows/*/from_bp": c_sup, "/rows/*/to_bp": c_sup,
                              "/rows/*/n_labels": c_sup + ": the set's size",
                              "/rows/*/n_labels_total": c_sup + ": the true count (no cut list on the arm)",
                              "/rows/*/joined": "DESIGN §6 who joined (" + c_sup + "); " + C_LABEL_LIST,
                              "/rows/*/left": "DESIGN §6 who left (" + c_sup + "); " + C_LABEL_LIST})
    clauses.update(ecl)
    if cut:
        cons.append("/rows/*/n_labels_total")
    partial = _cut_rows(ans, ctx, cons)
    red = [] if annotate else ["/total", "/rows/*/from_bp", "/rows/*/to_bp", "/rows/*/n_labels",
                               "/rows/*/n_labels_total", "/rows/*/joined", "/rows/*/left"]
    return _record(ans, clauses, cons + econ, partial=partial, rederived=red)


_TOOLS: Dict[str, Callable[[_Doc, dict, dict], dict]] = {
    "graphlet_summary": _graphlet_summary, "graphlet_walks": _graphlet_walks, "graphlet_walk": _graphlet_walk,
    "graphlet_labels": _graphlet_labels, "graphlet_splits": _graphlet_splits, "graphlet_sequence": _graphlet_sequence,
    "graphlet_export": _graphlet_export, "graphlet_subtrie": _graphlet_subtrie, "graphlet_save": _graphlet_save,
    "graphlet_load": _graphlet_load, "traverse_continue": _traverse_continue, "graphlet_claims": _graphlet_claims,
    "graphlet_support": _graphlet_support,
}


C_FAILED_SEED = "SPEC §7.1 (§7.0, §7.5): a seed whose permitted set could not be derived, or that a request budget " \
                "failed, has the failed shape {seed, outcome, error} with no arms and no graphlet; the code " \
                "no_graphlet is the runners' (FORMAT v2 B2)"


def _short(text: Any, n: int = 120) -> str:
    text = str(text)
    return text if len(text) <= n else text[:n] + "..."


def _failed_error(full_answer: dict, seed_index: int) -> Any:
    res = full_answer.get("results")
    r = res[seed_index] if isinstance(res, list) and seed_index < len(res) else full_answer
    return r.get("error")


def _arguments_refused(tool: str, args: dict, ctx: dict) -> Optional[dict]:
    """Arguments the tools refuse, or that W5 cannot judge, before any derivation:
    - `budget` on a server without local limits: bad_argument (library-defined, INVENTORY §1.1); with local limits a
      budget's stops are the tool's (DESIGN §5.4: deterministic, but library-defined): no claim (local_limits);
    (a `cursor` is judged in answer(), _cursor_record: its label selectors resolve against the document)."""
    if args.get("budget") is not None and tool not in NOT_APPLICABLE:
        if not ctx.get("local_limits"):
            return _record(None, {}, ["/*"], declared="expected_error", expected_error="bad_argument",
                           note="budget on a server without local limits (library-defined, INVENTORY §1.1)")
        return _declared("local_limits", "a budget under local limits: its stops are the tool's (DESIGN §5.4)")
    return None


# the documented defaults of the arguments a cursor binds (DESIGN §6's tool signatures; INVENTORY §2)
_DEFAULTS: Dict[str, dict] = {
    "graphlet_summary": {"arm": None, "detail": "summary", "budget": None},
    "graphlet_walks": {"rank": "support", "min_bp": 0, "label": None, "spell": "none", "tail_bp": 60,
                       "route_consistent": True, "coordinates": False, "budget": None},
    "graphlet_walk": {"budget": None},
    "graphlet_support": {"step": None, "budget": None},
    "graphlet_labels": {"arm": None, "rank": "direct_bp", "min_direct_bp": 0, "at_bp": None, "name": None,
                        "coordinates": False, "budget": None},
    "graphlet_splits": {"min_labels_before": 2, "budget": None},
    "graphlet_claims": {"min_bp": 0, "route_consistent": True, "labels": None, "budget": None},
    "graphlet_compare": {"arm": None, "mode": "claims", "budget": None},
}
_MISSING = object()


def _bound_args(d: "_Doc", tool: str, args: dict) -> Optional[dict]:
    """The arguments a cursor binds, normalised as the tool reads them (INVENTORY §1.2: all but cursor, max_bytes and
    n): arm short forms named (INVENTORY §1.5), label selectors resolved to labels (DESIGN §5 Graphlet.label: a
    label named by name, ref or id is one label), the documented defaults filled in (DESIGN §6). None where a
    selector does not resolve."""
    out = dict(_DEFAULTS.get(tool, {}))
    for k, v in args.items():
        if k in ("cursor", "max_bytes", "n", "access_token"):
            continue
        out[k] = v
    if isinstance(out.get("arm"), str):
        out["arm"] = _ARM_NAMES.get(out["arm"], out["arm"])
    try:
        for key in ("label", "name"):
            if out.get(key) is not None:
                out[key] = d.resolve(out[key])
        if out.get("labels") is not None:
            sel = out["labels"] if isinstance(out["labels"], list) else [out["labels"]]
            out["labels"] = sorted({d.resolve(x) for x in sel})
    except _Refused:
        return None
    return out


def _cursor_record(d: "_Doc", tool: str, args: dict, ctx: dict) -> dict:
    """A call that resumes a listing from a cursor (INVENTORY §1.2: opaque to W5, bound to the handle, the tool and
    the normalised arguments other than cursor, max_bytes and n). Its source (context cursor_source: the tool and
    arguments of the op that left it, found by the run step) decides: another tool, or an argument both name with
    another value after normalising (_bound_args), and the tool must refuse it: bad_cursor (library-defined,
    INVENTORY §1.2). Otherwise W5 makes no claim (cursor_page): a resumed page is the tool's, and an argument one
    call names without a documented default cannot be compared."""
    bad = _record(None, {}, ["/*"], declared="expected_error", expected_error="bad_cursor",
                  note="a cursor of another listing (library-defined, INVENTORY §1.2)")
    if ctx.get("cursor_bound_elsewhere"):
        return bad
    src = ctx.get("cursor_source")
    if not isinstance(src, dict) or not isinstance(src.get("args"), dict):
        return _declared("cursor_page", "a page resumed from an opaque cursor without a known source (INVENTORY §1.2)")
    if src.get("tool") != tool:
        return bad
    a, b = _bound_args(d, tool, src["args"]), _bound_args(d, tool, args)
    if a is None or b is None:
        return _declared("cursor_page", "a cursor whose label selectors do not resolve: the tool's")
    differs = False
    for k in set(a) | set(b):
        va, vb = a.get(k, _MISSING), b.get(k, _MISSING)
        if va is _MISSING or vb is _MISSING:
            return _declared("cursor_page", "%s: named by one call only, no documented default" % k)
        differs = differs or not same(va, vb)
    if differs:
        return bad
    return _declared("cursor_page", "a page resumed from the same listing: the tool's (INVENTORY §1.2)")


def _path_refused(args: dict, ctx: dict, spool: bool = False) -> Optional[dict]:
    """A file path of export, save or load (DESIGN §6: file names confined to the export directory; load also reads
    the spool, INVENTORY §2.17), judged on the path normalised and resolved against that directory:
    - inside, W5 answers (None): a plain file name once normalised (`sub/../a.fa` is `a.fa`), <export_dir>/<name>
      (<spool_dir>/<name> for load) as the runners write them, or an absolute path under the directory the context
      names (export_dir, spool_dir);
    - outside for sure, an expected error path_not_allowed (the code library-defined, INVENTORY §2.13, §2.17): the
      normalised path climbs out (`..`), or it is absolute and outside the directory (a runner writes any path inside
      the export or spool directory as <export_dir>/... or <spool_dir>/..., FORMAT v2 B1, so a raw absolute path is
      outside when the context names no directory);
    - no claim (path_unknown): `~` (the tool does not expand it: an io_error), the empty path, a subdirectory (whether
      it exists is the file system's)."""
    path = args.get("path")
    if path is None:
        return None
    unknown = _declared("path_unknown", "path %r: what the tool answers is its own" % (path,))
    if not isinstance(path, str) or not path or path.startswith("~"):
        return unknown
    out = _record(None, {}, ["/*"], declared="expected_error", expected_error="path_not_allowed",
                  note="path %r outside the export directory (DESIGN §6; the code library-defined, INVENTORY §2.13, "
                       "§2.17)" % path)
    symbols = ("<export_dir>", "<spool_dir>") if spool else ("<export_dir>",)
    roots = [r for r in (ctx.get("export_dir"), ctx.get("spool_dir") if spool else None) if isinstance(r, str)]

    def judge(rel: str) -> Optional[dict]:
        rel = posixpath.normpath(rel) if rel else "."
        if rel == ".." or rel.startswith("../"):
            return out
        return None if "/" not in rel and rel != "." else unknown

    for sym in symbols:
        if path == sym or path.startswith(sym + "/"):
            return judge(path[len(sym) + 1:])
    if path.startswith("<"):
        return unknown
    if posixpath.isabs(path):
        for root in roots:
            rel = posixpath.relpath(posixpath.normpath(path), posixpath.normpath(root))
            if not (rel == ".." or rel.startswith("../")):
                return judge(rel)
        return out
    return judge(path)


def seed_failed(full_answer: Optional[dict], seed_index: int = 0) -> bool:
    """Whether the full twin states the seed failed (SPEC §7.1: the result carries `error`, no arms)."""
    if not isinstance(full_answer, dict):
        return False
    res = full_answer.get("results")
    r = res[seed_index] if isinstance(res, list) and seed_index < len(res) else full_answer
    return isinstance(r, dict) and "error" in r


def answer(tool: str, args_resolved: dict, full_answer: dict, mgt_text: Optional[str],
           context: Optional[dict] = None) -> dict:
    """W5's record for one op (FORMAT-v8 §4); see the module docstring. On a seed the twin states failed, an op
    that reads the graphlet expects no_graphlet (there is none; |mgt_text| may be None)."""
    ctx = dict(context or {})
    args = {k: v for k, v in (args_resolved or {}).items() if k not in ("handle", "access_token", "cursor",
                                                                        "max_bytes", "budget")}
    if seed_failed(full_answer, ctx.get("seed_index", 0)) and (tool in _TOOLS or tool in PHASE_B or
                                                               tool == "graphlet_compare"):
        why = "the seed failed (%s): no graphlet" % _short(_failed_error(full_answer, ctx.get("seed_index", 0)))
        if tool == "graphlet_load":   # reads a file, not the graphlet: whatever the tool answers for no such file
            return _record(None, {}, ["/*"], declared="expected_error", expected_error=None,
                           note="%s (%s; the code is library-defined)" % (why, C_FAILED_SEED))
        return _record(None, {}, ["/*"], declared="expected_error", expected_error="no_graphlet",
                       note="%s (%s)" % (why, C_FAILED_SEED))
    if tool in PHASE_B:
        return _declared("not_implemented_phase_b", "%s: DESIGN-v8 §8 phase B" % tool)
    pre = _arguments_refused(tool, args_resolved or {}, ctx)
    if pre is not None:
        return pre
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
    if (args_resolved or {}).get("cursor") is not None:
        rec = _cursor_record(d, tool, args_resolved, ctx)
    else:
        try:
            rec = _TOOLS[tool](d, args, ctx)
        except _Refused as e:
            rec = _record(None, {}, ["/*"], declared="expected_error", expected_error=e.code,
                          note="%s (%s)" % (e, e.clause))
    gap = _memory_gap(d, ctx)
    if gap:   # every record read from this twin says so (VERIFY-D VD4)
        rec["oracle"]["inputs"] = list(rec["oracle"].get("inputs") or []) + [gap]
    return rec


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


COMPARABLE = (True, "qualified", "unverifiable", "unknown", False)   # DESIGN §5 compare's verdicts
_SPELLING_MODES = ("claims", "walks", "prefix_subset")


def compare_self_properties(ans: dict, full_answer: Optional[dict], mgt_text: str,
                            arm: Optional[str] = None) -> List[str]:
    """graphlet_compare(a, a) on the arms the compare covers (|arm|: the op's argument): the verdict rules on a
    graphlet compared with itself (DESIGN §5 compare, §5.2.5, §3.1; SPEC §7.5.7):
    - comparable is one of true, qualified, unverifiable, unknown, false; never false (one index identity, one
      oriented seed);
    - without index_fp in H: unverifiable ("unverifiable, never equal", DESIGN §3.1), or unknown only where a rule of
      unknown applies too: no bases in a mode that spells (DESIGN §5.2.5), or no common certified depth, complete_to_bp
      0 on a covered arm (library-defined: the tool's comparability, INVENTORY §2.14; DESIGN §5 only restricts a
      compare to min(complete_to_bp)), since DESIGN states no precedence between the two; with index_fp never
      unverifiable;
    - with index_fp, no bases in a mode that spells: unknown (DESIGN §5.2.5);
    - no difference under any verdict (counts zero or null, no rows); equal true exactly when comparable true; not
      true when the label evidence is a lower bound (DESIGN §5.2.5: what cannot be exact is qualified).
    -> problems."""
    d = _Doc(full_answer, mgt_text)
    sides = _compared_sides(d, d, arm)
    probs = []
    comp = ans.get("comparable")
    if not any(comp is v if isinstance(v, bool) else comp == v for v in COMPARABLE):
        probs.append("comparable %r: none of true, qualified, unverifiable, unknown, false (DESIGN §5)" % (comp,))
    if comp is False:
        probs.append("comparable false for a graphlet and itself (%s)" % ans.get("reason"))
    no_fp = d.h(14) is None
    spelling = ans.get("mode") in _SPELLING_MODES
    no_bases = spelling and any(d.flanks(s) is None for s in sides)
    depth0 = any(d.arms[s].get("complete_to_bp") == 0 for s in sides)
    if no_fp:
        allowed = ("unverifiable", "unknown") if (no_bases or depth0) else ("unverifiable",)
        if comp not in allowed or isinstance(comp, bool):
            probs.append("comparable %r without index_fp: %s (SPEC §7.5.7, DESIGN §3.1%s%s)" % (
                comp, " or ".join(allowed), ", §5.2.5" if no_bases else "",
                "; depth 0: library-defined, INVENTORY §2.14" if depth0 and not no_bases else ""))
    else:
        if comp == "unverifiable":
            probs.append("comparable unverifiable with index_fp set (SPEC §7.5.7, DESIGN §3.1)")
        if no_bases and comp != "unknown":
            probs.append("comparable %r without bases in mode %s (DESIGN §5.2.5: unknown)" % (comp, ans.get("mode")))
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


# ---------------------------------------------------------------------------------------------- the run step
W1_DIR = "dl-lib"          # FORMAT-v8 §3: W1's records, arms/dl-lib/graphlet_ops.jsonl
OUT_DIR = "oracle"         # W5's, arms/oracle/graphlet_ops.jsonl
_SYMBOL = re.compile(r"<handle(?::([^<>]+))?>\Z")


def _read_json(path: str) -> Any:
    with (gzip.open(path, "rt", encoding="utf-8") if path.endswith(".gz") else open(path, encoding="utf-8")) as fh:
        return json.load(fh)


def _inputs(run: str, sid: str, repeat: int) -> Tuple[Optional[dict], Optional[str], Optional[dict], dict]:
    """(the full twin, the MGT text, the graphlet response's envelope or None, the sources) of scenario |sid|: the
    twin at arms/direct-local/docs/<sid>.full.json.gz (FORMAT v2 C2); the engine's body from the direct-local raw
    answer of the graphlet request (VO12), else W1's document arms/dl-lib/docs/<sid>.mgt.gz (C1)."""
    safe = re.sub(r"[^A-Za-z0-9._-]", "_", str(sid))
    src: dict = {}
    full = text = envelope = None
    rel = "arms/direct-local/docs/%s.full.json.gz" % safe
    if os.path.isfile(os.path.join(run, rel)):
        full = _read_json(os.path.join(run, rel))
        full = full["answer"]["body"] if isinstance(full, dict) and isinstance(full.get("answer"), dict) else full
        src["full"] = rel
    rel = "arms/direct-local/" + suite_common.raw_name(sid, "direct-local", repeat)
    if os.path.isfile(os.path.join(run, rel)):
        obj = _read_json(os.path.join(run, rel))
        resp = obj["answer"]["body"] if isinstance(obj, dict) and isinstance(obj.get("answer"), dict) else obj
        res = (resp or {}).get("results") if isinstance(resp, dict) else None
        if isinstance(res, list) and res and isinstance(res[0], dict) and isinstance(res[0].get("graphlet"), str):
            text, envelope, src["mgt"] = res[0]["graphlet"], resp, rel + " (the engine's body)"
    if text is None:
        rel = "arms/%s/docs/%s.mgt.gz" % (W1_DIR, safe)
        if os.path.isfile(os.path.join(run, rel)):
            with gzip.open(os.path.join(run, rel), "rt", encoding="utf-8") as fh:
                text = fh.read()
            src["mgt"] = rel + " (W1's document)"
            # the graphlet request's normalised strategy, for _memory_gap and traverse_continue: the document's J
            # record is the response envelope (SPEC §7.5.2); W5 reads nothing else of it
            j = next((x[2:] for x in text.split("\n") if x.startswith("J ")), None)
            try:
                st = json.loads(j).get("strategy") if j else None
            except (ValueError, AttributeError):
                st = None
            if isinstance(st, dict):
                envelope, src["strategy"] = {"strategy": st}, rel + " J (the response envelope's strategy)"
    return full, text, envelope, src


def _context(rec: dict, view: bool, envelope: Optional[dict], full: dict) -> dict:
    """W5's context from W1's record (FORMAT-v8 §3 and v2 A4, B4): a listing W1 read to its end (paging_complete) is
    whole, else W5 checks W1's row count; the ceilings W1's library ran under; area G's one seed (seed_index 0)."""
    ans = rec.get("answer") if isinstance(rec.get("answer"), dict) else {}
    limits = (rec.get("library") or {}).get("limits") or {}
    args = rec.get("args_resolved") or {}
    ctx: dict = {"seed_index": 0, "paging": "all" if rec.get("paging_complete") is not False else "none",
                 "local_limits": limits.get("local_limits"),
                 "max_bytes": args.get("max_bytes", limits.get("max_bytes", 2048)),
                 "sequence_max_bytes": limits.get("sequence_max_bytes", 16384),
                 "strategy": (envelope or {}).get("strategy") or (full or {}).get("strategy")}
    if isinstance(ans.get("rows"), list):
        ctx["rows"] = len(ans["rows"])
    if view:
        ctx["view"] = True
    return ctx


_TORN_KEY = re.compile(r'"(scenario|op|tool)":\s*"((?:[^"\\]|\\.)*)"')


def _read_w1(path: str) -> Tuple[List[dict], List[dict], int]:
    """(records, torn, unkeyed) of W1's graphlet_ops.jsonl: the parsed lines; the torn lines whose scenario and op
    can still be read from the text (each becomes a W5 record declaring unreadable_w1); the number of torn lines
    with neither (one run-level record states them)."""
    recs, torn, unkeyed = [], [], 0
    if not os.path.exists(path):
        return recs, torn, unkeyed
    with open(path, encoding="utf-8", errors="replace") as fh:
        for no, line in enumerate(fh, 1):
            if not line.strip():
                continue
            try:
                obj = json.loads(line)
                if not isinstance(obj, dict):
                    raise ValueError("not an object")
                recs.append(obj)
            except ValueError:
                keys = dict(_TORN_KEY.findall(line))
                if keys.get("scenario") and keys.get("op"):
                    torn.append(dict(keys, line=no))
                else:
                    unkeyed += 1
    return recs, torn, unkeyed


def _target(hv: Any, sid: str, tools: Dict[Tuple[Any, Any], Any], scenarios: set) -> Tuple[Optional[str], bool]:
    """(the scenario whose graphlet a call reads, whether through a subtrie's view) from its handle symbol (FORMAT
    v2 B1): <handle> the scenario's, <handle:S> scenario S's, <handle:OP> / <handle:S/OP> the one op OP made (a
    subtrie's: a view; a load's: the same graphlet). (None, False) for a symbol naming no scenario or op (<handle:?>,
    a raw handle, a mistyped name)."""
    m = _SYMBOL.match(hv) if isinstance(hv, str) else None
    if not m:
        return None, False
    name = m.group(1)
    if name is None:
        return sid, False
    if "/" in name:
        owner, _, op = name.rpartition("/")
        return (owner, tools[(owner, op)] == "graphlet_subtrie") if (owner, op) in tools else (None, False)
    if (sid, name) in tools:
        return sid, tools[(sid, name)] == "graphlet_subtrie"
    return (name, False) if name in scenarios else (None, False)


def _cursor_source(rec: dict, w1: List[dict]) -> Optional[dict]:
    """The op a call's cursor comes from: the one earlier op of the scenario whose listing W1 did not read to its end
    (paging_complete false: the op that left a next_cursor for a selector); None when there is not exactly one."""
    earlier = []
    for r in w1:
        if r is rec:
            break
        if r.get("scenario") == rec.get("scenario") and r.get("paging_complete") is False and r.get("args_resolved"):
            earlier.append(r)
    return earlier[0] if len(earlier) == 1 else None


def oracle_records(run: str, repeat: int = 0, run_id: Optional[str] = None) -> List[dict]:
    """W5's records of a run (FORMAT v2 F4): one per W1 record of arms/dl-lib/graphlet_ops.jsonl, the same line
    shape (FORMAT-v8 §3) with way W5, W1's args_resolved, no library, and `oracle`; a torn W1 line gives a record
    declaring unreadable_w1 under its scenario and op where the text still names them, and the rest one run-level
    record (scenario and op null, run_level true) that the compare stage reports as inconclusive (R2-5)."""
    w1, torn, unkeyed = _read_w1(os.path.join(run, "arms", W1_DIR, "graphlet_ops.jsonl"))
    tools = {(r.get("scenario"), r.get("op")): r.get("tool") for r in w1}
    scenarios = {r.get("scenario") for r in w1}
    cache: Dict[str, tuple] = {}
    out = []
    for rec in w1:
        started = time.time()
        sid, tool = rec.get("scenario"), rec.get("tool")
        line = {"run_id": run_id or rec.get("run_id"), "scenario": sid, "way": "W5", "op": rec.get("op"),
                "tool": tool, "args_resolved": rec.get("args_resolved"), "selector_error": rec.get("selector_error"),
                "pages": None, "paging_complete": None, "answer": None, "error": None, "service_fields": None,
                "doc": None, "library": None, "elapsed_ms": 0.0, "started_utc": None, "sources": None}
        # whose graphlet the call reads (FORMAT v2 B1): <handle> the scenario's, <handle:S> scenario S's,
        # <handle:OP> / <handle:S/OP> the one an op made (a subtrie's view; a load's round trip of the same graphlet)
        hv = (rec.get("args_resolved") or {}).get("a" if tool == "graphlet_compare" else "handle", "<handle>")
        target, view = _target(hv, sid, tools, scenarios)
        if target is not None and target not in cache:
            cache[target] = _inputs(run, target, repeat)
        full, text, envelope, src = cache[target] if target is not None else (None, None, None, {})
        line["sources"] = src
        if rec.get("selector_error"):
            r = _declared("selector_error", "W1's arguments did not resolve: %s" % rec["selector_error"])
        elif rec.get("args_resolved") is None:
            r = _declared("no_input", "W1 recorded no arguments (%s)" % ((rec.get("error") or {}).get("code")))
        elif target is None:
            r = _declared("unresolved_handle", "the handle symbol %r names no scenario or op of this run" % (hv,))
        elif full is None or (text is None and not seed_failed(full)):
            r = _declared("no_input", "missing %s for scenario %s" % (
                " and ".join(k for k, v in (("the full twin", full), ("the MGT document", text)) if v is None), target))
        else:
            ctx = _context(rec, view, envelope, full)
            if (rec.get("args_resolved") or {}).get("cursor") is not None:
                src = _cursor_source(rec, w1)
                if src is not None:
                    ctx["cursor_source"] = {"tool": src.get("tool"), "args": src.get("args_resolved")}
            r = answer(tool, rec["args_resolved"], full, text, ctx)
        line["answer"], line["oracle"] = r["answer"], r["oracle"]
        line["elapsed_ms"] = round((time.time() - started) * 1000.0, 3)
        line["started_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(started))
        out.append(line)
    for t in torn:
        r = _declared("unreadable_w1", "W1's line %d is torn: its record cannot be read" % t["line"])
        out.append({"run_id": run_id, "scenario": t["scenario"], "way": "W5", "op": t["op"], "tool": t.get("tool"),
                    "args_resolved": None, "selector_error": None, "pages": None, "paging_complete": None,
                    "answer": None, "error": None, "service_fields": None, "doc": None, "library": None,
                    "elapsed_ms": 0.0, "started_utc": None, "sources": None, "oracle": r["oracle"]})
    if unkeyed:
        r = _declared("unreadable_w1", "%d torn W1 line(s) name no scenario and op: their ops have no W5 record"
                      % unkeyed)
        out.append({"run_id": run_id, "scenario": None, "way": "W5", "op": None, "tool": None, "run_level": True,
                    "args_resolved": None, "answer": None,
                    "error": {"code": "unreadable_w1_lines", "fields": {"lines": unkeyed}}, "oracle": r["oracle"]})
    return out


def write_run(run: str, repeat: int = 0, force: bool = False, quiet: bool = False) -> str:
    """Writes arms/oracle/graphlet_ops.jsonl (FORMAT v2 F4) and arms/oracle/run.json; an existing file is replaced
    only with |force| (never appended). -> the path written."""
    out_dir = os.path.join(run, "arms", OUT_DIR)
    path = os.path.join(out_dir, "graphlet_ops.jsonl")
    if os.path.exists(path) and not force:
        raise SystemExit("graphlet_oracle: %s exists (never appended): remove it or pass --force" % path)
    os.makedirs(out_dir, exist_ok=True)
    started = time.time()
    lines = oracle_records(run, repeat)
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        for line in lines:
            fh.write(json.dumps(line, sort_keys=True, ensure_ascii=False) + "\n")
    os.replace(tmp, path)
    counts: Dict[str, int] = {}
    for line in lines:
        k = (line.get("oracle") or {}).get("declared") or ("error" if line.get("error") else "answered")
        counts[k] = counts.get(k, 0) + 1
    summary = {"way": "W5", "records": len(lines), "outcomes": counts, "repeat": repeat,
               "seconds": round(time.time() - started, 3), "w1": "arms/%s/graphlet_ops.jsonl" % W1_DIR,
               "unreadable_w1_lines": sum(int(((x.get("error") or {}).get("fields") or {}).get("lines", 0))
                                          for x in lines if x.get("run_level"))
               + sum(1 for x in lines if (x.get("oracle") or {}).get("declared") == "unreadable_w1"
                     and not x.get("run_level"))}
    with open(os.path.join(out_dir, "run.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, sort_keys=True, indent=1)
    if not quiet:
        print("[W5] %d records %s -> %s" % (len(lines), json.dumps(counts, sort_keys=True), path))
    return path


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="W5: write arms/oracle/graphlet_ops.jsonl of a run directory from W1's "
                                             "records, the full twins and the graphlet documents (FORMAT v2 F4)")
    ap.add_argument("--run-dir", required=True, help="the run directory")
    ap.add_argument("--repeat", type=int, default=0, help="the direct-local repeat holding the graphlet answer")
    ap.add_argument("--force", action="store_true", help="replace an existing arms/oracle/graphlet_ops.jsonl")
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args(argv)
    write_run(a.run_dir, a.repeat, a.force, a.quiet)
    return 0


if __name__ == "__main__":
    sys.exit(main())
