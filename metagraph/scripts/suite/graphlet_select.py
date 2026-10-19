"""Run-time selectors of area G (DESIGN-v8 §4, owner decision D4; FORMAT-v8 §2): an op's argument named by what an
earlier op of the same way answered, resolved in each way from that way's own answers, so that no argument is frozen
before the run and a selector that resolves differently in two ways is named by the compare stage (args_resolved).

    resolve(value, answers, where="args") -> the literal value to send
        value:   a literal, {"ref": "<op id>", "path": "<RFC 6901 pointer>"}, or {"select": "<name>", "args": {...}},
                 at any depth inside lists and objects (a literal object is walked; only these exact key sets are
                 selectors: a label selector {"ref": "c:3"} has no "path" and stays a literal)
        answers: {op id: the op's ASSEMBLED answer in this way (FORMAT-v8 §3: list rows concatenated over its pages),
                 or None when the op failed, its own selector failed, or it does not run in this way}
    SELECTORS: {name: function(answers, **args)}; each defined in the suite's words in its docstring

A selector that cannot answer raises SelectorError with a message naming the argument, the op and the rule; the
runner records it as the op's `selector_error` and makes no call. Pure: no I/O, no library import, stdlib only,
Python 3.8+; async's MCP runner vendors this file unchanged (as it does checks.py).

The answers a selector reads are the library's tool answers as the suite sees them (metagraph/api/python/metagraph/
traverse/mcp_tools.py): every paged list's rows are `rows`; a walk is named by its path id `walk` (graphlet_walks rows,
graphlet_walk's base); a label is {name, ref}, `ref` its cross-call key. A selector returns a label as the label
selector {"ref": ref} by default (form="ref" the ref string, form="name" the name), the form every label argument of
the tools takes.
"""

from __future__ import annotations

import copy
import inspect
import re
from typing import Any, Callable, Dict, List, Mapping, Optional

ARMS = ("left", "right")
LABEL_FORMS = ("selector", "ref", "name")
_INDEX = re.compile(r"(0|[1-9][0-9]*)\Z")


class SelectorError(Exception):
    """A selector that cannot answer from the earlier answers of its way: the op is not sent (selector_error)."""


# ------------------------------------------------------------------------------------------------- resolution
def is_ref(value: Any) -> bool:
    """{"ref": "<op id>", "path": "<pointer>"}: exactly these two keys."""
    return isinstance(value, dict) and set(value) == {"ref", "path"}


def is_select(value: Any) -> bool:
    """{"select": "<name>"} or {"select": "<name>", "args": {...}}: exactly these keys."""
    return isinstance(value, dict) and set(value) in ({"select"}, {"select", "args"})


def resolve(value: Any, answers: Mapping[str, Any], where: str = "args") -> Any:
    """|value| with every op reference and named selector replaced by what it selects from |answers| (a new object;
    |value| and |answers| are not changed). |where| names the argument in an error's message."""
    if is_ref(value):
        op, path = value["ref"], value["path"]
        if not isinstance(op, str) or not isinstance(path, str):
            raise SelectorError("%s: a reference is {\"ref\": \"<op id>\", \"path\": \"<JSON pointer>\"}, not %r"
                                % (where, value))
        try:
            return copy.deepcopy(pointer(answer_of(answers, op), path))
        except SelectorError as e:
            raise SelectorError("%s: %s" % (where, e)) from None
    if is_select(value):
        name, args = value["select"], value.get("args", {})
        if name not in SELECTORS:
            raise SelectorError("%s: no named selector %r (known: %s)" % (where, name, ", ".join(sorted(SELECTORS))))
        if not isinstance(args, dict):
            raise SelectorError("%s: the args of selector %s are an object, not %r" % (where, name, args))
        args = dict((k, resolve(v, answers, "%s.%s.args.%s" % (where, name, k))) for k, v in args.items())
        fn = SELECTORS[name]
        try:
            inspect.signature(fn).bind(answers, **args)
        except TypeError as e:   # an argument the selector does not take, or one it needs left out
            raise SelectorError("%s: %s: %s" % (where, name, e)) from None
        try:
            out = fn(answers, **args)
        except SelectorError as e:
            raise SelectorError("%s: %s: %s" % (where, name, e)) from None
        return copy.deepcopy(out)
    if isinstance(value, dict):
        return dict((k, resolve(v, answers, "%s.%s" % (where, k))) for k, v in value.items())
    if isinstance(value, list):
        return [resolve(v, answers, "%s[%d]" % (where, i)) for i, v in enumerate(value)]
    return value


def answer_of(answers: Mapping[str, Any], op: Any) -> Any:
    """The assembled answer of earlier op |op| in this way; SelectorError when there is none to read."""
    if not isinstance(op, str):
        raise SelectorError("an op is named by its id (a string), not %r" % (op,))
    if op not in answers:
        raise SelectorError("no op %r before this one in this way" % op)
    ans = answers[op]
    if ans is None:
        raise SelectorError("op %r has no answer in this way (it failed, its own selector failed, or it does not "
                            "run in this way)" % op)
    if isinstance(ans, dict) and isinstance(ans.get("error"), str):
        raise SelectorError("op %r answered the error %s" % (op, ans["error"]))
    return ans


def pointer(doc: Any, path: str) -> Any:
    """RFC 6901: the value at JSON pointer |path| in |doc| ("" the whole document; "/a/0/b"; "~1" is "/", "~0" is
    "~"). A path that does not lead to a value is a SelectorError naming the step."""
    if path == "":
        return doc
    if not path.startswith("/"):
        raise SelectorError("the JSON pointer %r does not start with '/'" % path)
    cur = doc
    seen = ""
    for raw in path[1:].split("/"):
        if re.search(r"~(?![01])", raw):
            raise SelectorError("the JSON pointer %r: '~' is escaped as ~0 or ~1 only" % path)
        tok = raw.replace("~1", "/").replace("~0", "~")
        if isinstance(cur, dict):
            if tok not in cur:
                raise SelectorError("the JSON pointer %r: no key %r at %r" % (path, tok, seen or "/"))
            cur = cur[tok]
        elif isinstance(cur, list):
            if not _INDEX.match(tok):
                raise SelectorError("the JSON pointer %r: %r at %r is not an index of a list" % (path, tok, seen or "/"))
            if int(tok) >= len(cur):
                raise SelectorError("the JSON pointer %r: index %s at %r, the list has %d items"
                                    % (path, tok, seen or "/", len(cur)))
            cur = cur[int(tok)]
        else:
            raise SelectorError("the JSON pointer %r: %r at %r is a %s, not an object or a list"
                                % (path, tok, seen or "/", type(cur).__name__))
        seen += "/" + raw
    return cur


# --------------------------------------------------------------------------------------------- answer readers
def _ops(names: Any, what: str) -> List[str]:
    """An op id or a list of op ids (at least one)."""
    if isinstance(names, str):
        return [names]
    if isinstance(names, list) and names and all(isinstance(x, str) for x in names):
        return list(names)
    raise SelectorError("%s names an op id or a non-empty list of op ids, not %r" % (what, names))


def _walk_answer(answers: Mapping[str, Any], op: str) -> dict:
    """A graphlet_walk answer: {arm, walk, length_bp, rows: [{segment, in, out, labels_out, merge_of?, ...}],
    continuation?} (mcp_tools.py graphlet_walk, _walk_base, _walk_rows)."""
    ans = answer_of(answers, op)
    if not (isinstance(ans, dict) and isinstance(ans.get("walk"), int) and isinstance(ans.get("rows"), list)
            and all(isinstance(r, dict) and "segment" in r for r in ans["rows"])):
        raise SelectorError("op %r is not a graphlet_walk answer (walk and its segment rows)" % op)
    return ans


def _int_arg(name: str, value: Any, least: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < least:
        raise SelectorError("%s is an integer >= %d, not %r" % (name, least, value))
    return value


def _choice(name: str, value: Any, allowed: tuple) -> Any:
    if value not in allowed:
        raise SelectorError("%s is one of %s, not %r" % (name, ", ".join(repr(x) for x in allowed), value))
    return value


def _label(row: dict, form: str) -> Any:
    _choice("form", form, LABEL_FORMS)
    if not isinstance(row, dict) or not isinstance(row.get("ref"), str):
        raise SelectorError("a label is {name, ref}, not %r" % (row,))
    if form == "ref":
        return row["ref"]
    if form == "name":
        return row.get("name")
    return {"ref": row["ref"]}


def _pick_walk(answers: Mapping[str, Any], walks: Any, test: Callable[[dict], bool], what: str, order: str,
               arm: Optional[str]) -> int:
    """The walk id of the graphlet_walk answers named by |walks| that pass |test|: the lowest walk id (order "id")
    or the first in the order named (order "given"), of arm |arm| when given."""
    _choice("order", order, ("id", "given"))
    if arm is not None:
        _choice("arm", arm, ARMS)
    names = _ops(walks, "walks")
    hits = []
    for op in names:
        ans = _walk_answer(answers, op)
        if arm is not None and ans.get("arm") != arm:
            continue
        if test(ans):
            hits.append(ans)
    if not hits:
        raise SelectorError("none of the walks of %s%s %s" % (", ".join(names), "" if arm is None else
                                                               " on the %s arm" % arm, what))
    if arm is None and len(set(h.get("arm") for h in hits)) > 1:
        raise SelectorError("the walks that %s lie on both arms: name the arm" % what)
    if order == "id":
        hits.sort(key=lambda h: h["walk"])
    return hits[0]["walk"]


def _merge_rows(ans: dict) -> List[int]:
    """The chain positions of a graphlet_walk answer's merge segments (rows with merge_of: two or more parents)."""
    return [i for i, r in enumerate(ans["rows"]) if isinstance(r.get("merge_of"), list) and len(r["merge_of"]) > 1]


# --------------------------------------------------------------------------------------------- named selectors
def walk_through_merge(answers: Mapping[str, Any], walks: Any, order: str = "id", arm: Optional[str] = None) -> int:
    """A walk through a merge: of the walks whose graphlet_walk answers |walks| names (an op id or a list of them,
    each assembled whole: paging "all"), one whose segment chain holds a merge, i.e. a row with `merge_of` listing
    two or more parents: the chain is the walk's first-parent path from its leaf to the root, so the walk reaches
    that segment from both parents' branches. The lowest such walk id (order "id", the default) or the first in the
    order named (order "given"), on arm |arm| when given. -> the walk id."""
    return _pick_walk(answers, walks, lambda a: bool(_merge_rows(a)), "passes a merge (a row with merge_of)",
                      order, arm)


def walk_with_continuation(answers: Mapping[str, Any], walks: Any, order: str = "id", arm: Optional[str] = None,
                           min_labels: int = 0) -> int:
    """A walk with a continuation: of the walks whose graphlet_walk answers |walks| names, one whose answer carries
    a `continuation` block (its leaf records the bases and labels a traverse_continue of that walk starts from)
    with at least |min_labels| labels. The lowest such walk id (order "id") or the first named (order "given"), on
    arm |arm| when given. -> the walk id."""
    least = _int_arg("min_labels", min_labels, 0)

    def has(a: dict) -> bool:
        c = a.get("continuation")
        return isinstance(c, dict) and isinstance(c.get("labels"), int) and c["labels"] >= least
    return _pick_walk(answers, walks, has, "has a continuation%s" % (" with %d labels or more" % least
                                                                       if least else ""), order, arm)


def slice_past_end(answers: Mapping[str, Any], length_of: str, part: str, walk: Optional[int] = None,
                   keep: int = 0, past: int = 1) -> int:
    """One end of a slice [start, end) of a walk's bases that runs past the walk's last base: graphlet_sequence
    clips such an end to the bases there are (the answer's `to`), which this slice tests. L is the walk's length as
    op |length_of| answered it: a graphlet_sequence answer's `length` (the bases that call spelled, the seed with
    with_seed), a graphlet_walk answer's `length_bp`, or the `length_bp` of the graphlet_walks row whose walk is
    |walk| (required for a list). part "start" -> max(0, L - keep) (keep 0: the slice starts at the end, and is
    empty); part "end" -> L + past (past >= 1: beyond the last base)."""
    _choice("part", part, ("start", "end"))
    keep = _int_arg("keep", keep, 0)
    past = _int_arg("past", past, 1)
    ans = answer_of(answers, length_of)
    if not isinstance(ans, dict):
        raise SelectorError("op %r is not a graphlet_sequence, graphlet_walk or graphlet_walks answer" % length_of)
    if isinstance(ans.get("sequence"), str) and isinstance(ans.get("length"), int):
        length = ans["length"]
    elif isinstance(ans.get("walk"), int) and isinstance(ans.get("length_bp"), int):
        if walk is not None and walk != ans["walk"]:
            raise SelectorError("op %r answered walk %d, not walk %r" % (length_of, ans["walk"], walk))
        length = ans["length_bp"]
    elif isinstance(ans.get("rows"), list):
        if walk is None:
            raise SelectorError("op %r is a list of walks: name the walk" % length_of)
        rows = [r for r in ans["rows"] if isinstance(r, dict) and r.get("walk") == walk]
        if not rows or not isinstance(rows[0].get("length_bp"), int):
            raise SelectorError("op %r lists no walk %r with its length_bp" % (length_of, walk))
        length = rows[0]["length_bp"]
    else:
        raise SelectorError("op %r is not a graphlet_sequence, graphlet_walk or graphlet_walks answer" % length_of)
    return max(0, length - keep) if part == "start" else length + past


def label_entered_at_merge(answers: Mapping[str, Any], walk: Any, form: str = "selector") -> Any:
    """A label that entered a walk at a merge: in the graphlet_walk answer(s) |walk| names (the first named first),
    at the first merge row of the chain (a row with merge_of), the first label of the row's `in` list, in the order
    shown, that the previous row's `out` does not hold. The previous row is the merge's first parent (the chain
    follows first parents), so such a label is alive at the merge only because another parent carried it in: the
    displayed walk shows it from the merge on (a merge-entered label). The previous row's `out` must be shown whole
    (its list as long as labels_out), else absence there proves nothing and the selector fails; a merge whose shown
    `in` labels all come from the first parent is passed over. -> the label (form "selector" {"ref"}, "ref",
    "name")."""
    _choice("form", form, LABEL_FORMS)
    names = _ops(walk, "walk")
    seen = []
    for op in names:
        ans = _walk_answer(answers, op)
        rows = ans["rows"]
        for i in _merge_rows(ans):
            row = rows[i]
            if i == 0:
                continue   # a merge at the chain's start has no first-parent row in this answer
            prev = rows[i - 1]
            if prev.get("segment") != row["merge_of"][0]:
                raise SelectorError("op %r: the row before merge segment %r is segment %r, not its first parent %r"
                                    % (op, row.get("segment"), prev.get("segment"), row["merge_of"][0]))
            out = prev.get("out")
            if not isinstance(out, list) or not isinstance(prev.get("labels_out"), int):
                raise SelectorError("op %r: segment %r's row has no out list with labels_out"
                                    % (op, prev.get("segment")))
            if len(out) < prev["labels_out"]:
                raise SelectorError("op %r: segment %r shows %d of its %d out labels: a label missing there may "
                                    "still be one of them" % (op, prev.get("segment"), len(out), prev["labels_out"]))
            have = set(l.get("ref") for l in out if isinstance(l, dict))
            for lab in row.get("in") or []:
                if isinstance(lab, dict) and lab.get("ref") not in have:
                    return _label(lab, form)
            seen.append("%s segment %s" % (op, row.get("segment")))
    if seen:
        raise SelectorError("every label shown entering %s also leaves the first parent" % ", ".join(seen))
    raise SelectorError("no merge on the walks of %s" % ", ".join(names))


def top_label(answers: Mapping[str, Any], labels: str, arm: Optional[str] = None, index: int = 0,
              form: str = "selector") -> Any:
    """The top label of a label table: row |index| (0: the first) of the graphlet_labels table answer |labels|
    (rows {name, ref, arm, direct_bp, reach_bp, ...}, ranked by the call's rank, then ref, then arm), counting only the
    rows of arm |arm| when given. -> the label (form "selector" {"ref"}, "ref", "name")."""
    index = _int_arg("index", index, 0)
    if arm is not None:
        _choice("arm", arm, ARMS)
    ans = answer_of(answers, labels)
    rows = ans.get("rows") if isinstance(ans, dict) else None
    if not isinstance(rows, list) or "name" in ans or \
            not all(isinstance(r, dict) and "direct_bp" in r and "ref" in r for r in rows):
        raise SelectorError("op %r is not a graphlet_labels table answer (rows with ref and direct_bp)" % labels)
    rows = [r for r in rows if arm is None or r.get("arm") == arm]
    if index >= len(rows):
        raise SelectorError("op %r lists %d label rows%s, no row %d" % (labels, len(rows), "" if arm is None else
                                                                          " on the %s arm" % arm, index))
    return _label(rows[index], form)


def _arms_of(answers: Mapping[str, Any], of: str) -> List[str]:
    """The arms a graphlet answers for, in the order left, right: the keys of a graphlet_summary(detail="summary",
    arm=None) answer's summary.arms, else those of its evidence block's complete_to_bp."""
    ans = answer_of(answers, of)
    arms = None
    if isinstance(ans, dict) and isinstance(ans.get("summary"), dict) and isinstance(ans["summary"].get("arms"),
                                                                                       dict):
        arms = ans["summary"]["arms"]
    elif isinstance(ans, dict) and isinstance(ans.get("evidence"), dict) and \
            isinstance(ans["evidence"].get("complete_to_bp"), dict):
        arms = ans["evidence"]["complete_to_bp"]
    if not arms or any(a not in ARMS for a in arms):
        raise SelectorError("op %r names no arms (a graphlet_summary answer's summary.arms, or an evidence block's "
                            "complete_to_bp)" % of)
    return [a for a in ARMS if a in arms]


def first_arm(answers: Mapping[str, Any], of: str) -> str:
    """The graphlet's first arm in the order left, right, as the graphlet_summary(detail="summary", arm=None) answer
    |of| lists its arms: "left" when the graphlet has a left arm, else "right" (a right-only graphlet). With
    second_arm_or_none, an op plan names both arms of any graphlet without knowing which it has. -> "left" or
    "right"."""
    return _arms_of(answers, of)[0]


def second_arm_or_none(answers: Mapping[str, Any], of: str) -> Optional[str]:
    """The graphlet's second arm, or none: "right" when the graphlet_summary(detail="summary", arm=None) answer |of|
    lists both arms, None when it lists one (a one-arm graphlet, where every tool takes arm=None for its only arm).
    -> "right" or None."""
    arms = _arms_of(answers, of)
    return arms[1] if len(arms) > 1 else None


SELECTORS: Dict[str, Callable[..., Any]] = {
    "walk_through_merge": walk_through_merge,
    "walk_with_continuation": walk_with_continuation,
    "slice_past_end": slice_past_end,
    "label_entered_at_merge": label_entered_at_merge,
    "top_label": top_label,
    "first_arm": first_arm,
    "second_arm_or_none": second_arm_or_none,
}
