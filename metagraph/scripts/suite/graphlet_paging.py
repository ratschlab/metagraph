"""Paging and assembly of area G (FORMAT-v8 §1, §3; FORMAT v2 A1-A4): how a runner follows a graphlet tool's pages and
makes ONE answer of them, shared by every runner so that two ways never assemble differently. Lane L's local runner
(run_graphlet_local.py) imports it; async's MCP runner vendors it unchanged (as it does graphlet_select.py). Pure: no
I/O, no library import, stdlib only, Python 3.8+.

    paging_problem(paging)            -> None, or why an op's `paging` is malformed ("none" | "all" | {"pages": N >= 1})
    max_pages(paging)                 -> the pages a paging mode fetches at most
    error_page(page)                  -> whether a page is a tool error ({"error": "<code>", ...})
    error_record(page)                -> FORMAT's {"code", "fields"} of an error page (its text left out)
    next_call(tool, page)             -> (argument, value) of the next page's call, or None after the last page
    more(tool, paging, pages, given)  -> whether the loop fetches another page
    loop_error(tool, paging, pages, given) -> None, or the paging_runaway record of a loop that must stop
    assemble(tool, pages, complete, volatile=())
                                      -> the pages as one answer (AssemblyError when they do not make one)
    paging_complete(tool, last, answer) -> FORMAT v2 A4's record field: whether the listing was followed to its end
    selector_view(tool, answer, last) -> the answer later selectors of the way read (graphlet_select's `answers`)

The library's tools (metagraph/api/python/metagraph/traverse/mcp_tools.py) page every list as `rows` with a
`next_cursor` (GraphletTools._page :989-1085, the page at :1008-1013; graphlet_compare by _page_cut :1087-1123, whose
base fields give way per page, named in fields_cut), and graphlet_sequence by `next_from` (:2150-2161: a slice longer
than the ceiling is cut; the next page is the same call with start = next_from). A runner calls, for one op:

    given = the op's own continuation argument, ("cursor", args["cursor"]) or ("start", args["start"]), or None
    page = call(args); pages = [page]
    while not error_page(page) and more(tool, paging, pages, given):
        arg, value = next_call(tool, page); page = call(dict(args, **{arg: value})); pages.append(page)
    error page: record error_record(page); else loop_error(tool, paging, pages, given): record that error;
    else answer = assemble(tool, pages, complete=next_call(tool, pages[-1]) is None)
         record["paging_complete"] = paging_complete(tool, pages[-1], answer)
         answers[op id] = selector_view(tool, record answer, pages[-1])   # (None for an op with an error)

(an op with an error has no answer and paging_complete null).
"""

from __future__ import annotations

from typing import Any, List, Optional, Tuple

# A1: the list field of each paged tool; tests/test_graphlet_paging.py checks the table against the library's source
# (every self._page / self._page_cut call names one of these tools)
LIST_FIELDS = {
    "graphlet_summary": "rows",    # detail="limitations" (:1595); detail="summary" is one object
    "graphlet_walks": "rows",      # :1676
    "graphlet_walk": "rows",       # :1802; budgeted :1748
    "graphlet_support": "rows",    # :1838
    "graphlet_labels": "rows",     # table :1930, name= :1896; budgeted name= :1961
    "graphlet_splits": "rows",     # :1986; budgeted :2018
    "graphlet_claims": "rows",     # :2094; budgeted :2097
    "graphlet_compare": "rows",    # _page_cut :2361
    "graphlet_list": "rows",       # :2400
}
SEQUENCE_TOOL = "graphlet_sequence"
# the fields that say how a listing ended (a budgeted list's total_at_least, complete, filtered_so_far, its local block
# and stop, :1008-1014; the hint that goes with `filtered`, walks :1669-1675, claims :2087-2089, which a budgeted
# claims listing states only once its resumed pages have filtered something): the last page's
END_FIELDS = ("total", "total_at_least", "complete", "filtered", "filtered_so_far", "hint", "local", "stop")
# fields of one page's fitting, merged over the pages; `cursor` (FORMAT §3: "next_cursor/cursor removed") is a page's
# own echo where a service states one
PAGE_FIELDS = ("next_cursor", "cursor", "fields_cut", "row_truncated", "cut_fields")
# A3: the marker on each row of a page whose one row was cut to fit the page alone (_page :1068-1079: row_truncated
# and cut_fields; _shrink :2713-2733: a string loses its tail, a list its second half, scalars never): that page's own
# cut_fields, so that a reader knows which assembled row is a cut of the uncut row
ROW_CUT_MARKER = "_cut_fields"
# an error's text, kept in the raw page only (FORMAT v2 B2: errors by code and structured fields)
ERROR_TEXT = ("error", "message", "hint")
# graphlet_select.CURSOR_KEY (the same string): where selector_view keeps a stopped listing's last next_cursor
SELECTOR_CURSOR_KEY = "_next_cursor"
# "all" stops here: a listing that never ends is a defect, recorded as paging_runaway (loop_error)
MAX_PAGES = 10000


class AssemblyError(Exception):
    """Pages that do not make one answer (a list page without its list, base fields two pages state differently)."""


def paging_problem(paging: Any) -> Optional[str]:
    """None for a valid op paging ("none": one call; "all": every page; {"pages": N}: at most N, N >= 1), else why."""
    if paging in ("none", "all"):
        return None
    if isinstance(paging, dict) and set(paging) == {"pages"} and isinstance(paging["pages"], int) \
            and not isinstance(paging["pages"], bool) and paging["pages"] >= 1:
        return None
    return "paging is \"none\", \"all\" or {\"pages\": N >= 1}, not %r" % (paging,)


def max_pages(paging: Any) -> int:
    if paging == "none":
        return 1
    if paging == "all":
        return MAX_PAGES
    return paging["pages"]


def error_page(page: Any) -> bool:
    return isinstance(page, dict) and isinstance(page.get("error"), str)


def error_record(page: dict) -> dict:
    """An error page -> {"code", "fields"}: its structured fields (replayable, reason, stop, of, max_bytes, ...), never
    its message or hint text."""
    return {"code": page["error"], "fields": dict((k, v) for k, v in page.items() if k not in ERROR_TEXT)}


def next_call(tool: str, page: Any) -> Optional[Tuple[str, Any]]:
    """(the argument, its value) that asks for the page after |page|, or None when |page| is the last: graphlet_sequence
    continues with start = next_from, every other tool with cursor = next_cursor."""
    if not isinstance(page, dict):
        return None
    if tool == SEQUENCE_TOOL:
        return None if page.get("next_from") is None else ("start", page["next_from"])
    return None if page.get("next_cursor") is None else ("cursor", page["next_cursor"])


def _repeats(tool: str, pages: List[Any], given: Optional[Tuple[str, Any]]) -> bool:
    """Whether the last page's continuation asks for a page already asked for: the cursor or next_from that asked for
    the last page (|given| for the first), and for graphlet_sequence a next_from not beyond the page's own `from`."""
    nxt = next_call(tool, pages[-1])
    if nxt is None:
        return False
    before = next_call(tool, pages[-2]) if len(pages) > 1 else given
    if before is not None and before == nxt:
        return True
    start = pages[-1].get("from") if isinstance(pages[-1], dict) else None
    return tool == SEQUENCE_TOOL and isinstance(start, int) and isinstance(nxt[1], int) and nxt[1] <= start


def loop_error(tool: str, paging: Any, pages: List[Any], given: Optional[Tuple[str, Any]] = None) -> Optional[dict]:
    """The paging_runaway record of a loop that must stop although its paging asks for more (VERIFY-L VL6), else None:
    the last page names a next page that repeats one already asked for (a cursor or next_from that does not advance:
    the same page again, forever), or paging "all" fetched MAX_PAGES pages and a next one is still named. A `pages`
    cap reached, or paging "none", is no runaway: the op asked for no more."""
    if not pages or next_call(tool, pages[-1]) is None or len(pages) >= max_pages(paging) and paging != "all":
        return None
    if _repeats(tool, pages, given):
        return {"code": "paging_runaway", "fields": {"pages": len(pages), "reason": "repeats"}}
    if paging == "all" and len(pages) >= MAX_PAGES:
        return {"code": "paging_runaway", "fields": {"pages": len(pages), "reason": "max_pages"}}
    return None


def more(tool: str, paging: Any, pages: List[Any], given: Optional[Tuple[str, Any]] = None) -> bool:
    """Whether the loop fetches another page: the last page names one, the paging allows it, and it is no runaway."""
    return next_call(tool, pages[-1]) is not None and len(pages) < max_pages(paging) \
        and loop_error(tool, paging, pages, given) is None


def unmarked(evidence: Any) -> Any:
    """An evidence block without the local-stop mark a page its budget stopped carries (limitations.local,
    mcp_tools._mark_interrupted :248-262); the page that completes the listing has none."""
    if isinstance(evidence, dict) and isinstance(evidence.get("limitations"), dict) \
            and "local" in evidence["limitations"]:
        lims = dict(evidence["limitations"])
        lims.pop("local")
        return dict(evidence, limitations=lims)
    return evidence


def assemble(tool: str, pages: List[Any], complete: bool, volatile: Tuple[str, ...] = ()) -> Any:
    """FORMAT v2 A2: the pages of one op as one answer. A list tool's `rows` concatenated in page order, every other
    field from the first page that states it (page 0's, and a key page 0 lacks from the first later page carrying it:
    nothing a page states is dropped), except: the END_FIELDS from the last page; `evidence` from the last page that
    states it (an interrupted page's carries limitations.local, the completing page's not); compare's base fields that
    a page cut (fields_cut) thereby taken from any page that carries them, fields_cut then naming only those no page
    carried;
    `row_truncated` / `cut_fields` (a row cut to fit its page alone) kept when any page has them, the fields merged,
    and each row of such a page marked "_cut_fields": that page's cut_fields (A3); next_cursor and cursor removed.
    |volatile|: keys a runner's pages may state differently page by page, taken from the last page that states them
    and exempt from the check below. The library's own pages have none (the local runner passes none); a service
    that adds per-page fields of its own (an elapsed time, a lane, an op id) names them here, and only them, so that
    a field the library states differently on two pages is still an AssemblyError (VERIFY-L VL4).
    graphlet_sequence: `sequence` concatenated over the next_from loop, `to` from the last page,
    next_from removed and `truncated` too when the loop reached the end (|complete|). A base field two pages state
    differently (the evidence but for that mark) is an AssemblyError: a listing's pages describe one answer. The pages
    are not changed."""
    if not pages:
        raise AssemblyError("no page of %s to assemble" % tool)
    first = dict(pages[0]) if isinstance(pages[0], dict) else pages[0]
    if not isinstance(first, dict):
        if len(pages) > 1:
            raise AssemblyError("%d pages of %s that are no objects" % (len(pages), tool))
        return first
    first.pop("next_cursor", None)
    first.pop("cursor", None)
    if tool == SEQUENCE_TOOL:
        if not all(isinstance(p, dict) and isinstance(p.get("sequence"), str) for p in pages):
            raise AssemblyError("a page of graphlet_sequence without its sequence")
        first["sequence"] = "".join(p["sequence"] for p in pages)
        first["to"] = pages[-1].get("to")
        first.pop("next_from", None)
        if complete:
            first.pop("truncated", None)
        return first
    field = LIST_FIELDS.get(tool)
    if field is None or field not in first:
        if len(pages) > 1:
            raise AssemblyError("%d pages of %s, the first without a list" % (len(pages), tool))
        return first
    rows = []
    for i, p in enumerate(pages):
        if not isinstance(p, dict) or not isinstance(p.get(field), list):
            raise AssemblyError("page %d of %s has no %s list" % (i, tool, field))
        if p.get("row_truncated"):
            mark = list(p.get("cut_fields") or [])
            rows.extend(dict(r, **{ROW_CUT_MARKER: mark}) if isinstance(r, dict) else r for r in p[field])
        else:
            rows.extend(p[field])
        for k, v in p.items():
            if k in END_FIELDS or k in PAGE_FIELDS or k == field or k in volatile:
                continue
            if k not in first:   # page 0 lacks it (cut there, or stated only later): the first page stating it
                first[k] = v
            elif unmarked(v) != unmarked(first[k]) if k == "evidence" else v != first[k]:
                raise AssemblyError("page %d of %s states %s otherwise than an earlier page" % (i, tool, k))
    first[field] = rows
    for k in ("evidence",) + tuple(volatile):
        for p in reversed(pages):
            if k in p:
                first[k] = p[k]
                break
    last = pages[-1]
    for k in END_FIELDS:
        if k in last:
            first[k] = last[k]
        else:
            first.pop(k, None)
    cut = []
    for p in pages:
        for k in p.get("fields_cut") or []:
            if k not in cut and k not in first:
                cut.append(k)
    if cut:
        first["fields_cut"] = cut
    else:
        first.pop("fields_cut", None)
    cut_fields = []
    for p in pages:
        if p.get("row_truncated"):
            cut_fields += [k for k in p.get("cut_fields") or [] if k not in cut_fields]
    if cut_fields or any(p.get("row_truncated") for p in pages):
        first["row_truncated"] = True
        first["cut_fields"] = cut_fields
    else:
        first.pop("row_truncated", None)
        first.pop("cut_fields", None)
    return first


def paging_complete(tool: str, last: Any, answer: Any) -> bool:
    """FORMAT v2 A4: whether the op followed its listing to the end: false when the last page fetched still names a
    next page (a `pages` cap, paging "none" on a longer list) or the assembled answer states complete false (a stop
    of its local budget the pages did not get past); true for an answer that is no listing."""
    return next_call(tool, last) is None and not (isinstance(answer, dict) and answer.get("complete") is False)


def selector_view(tool: str, answer: Any, last: Any) -> Any:
    """The answer the way's later selectors read for an op (graphlet_select.resolve's `answers`): the record's answer,
    plus, when the op's listing stopped before its end with a cursor left (paging "none" on a longer list, a `pages`
    cap), that last page's next_cursor under SELECTOR_CURSOR_KEY, which graphlet_select's {"cursor_of": op} reads
    (VERIFY-L VL5). The record's answer itself never carries it (FORMAT §3)."""
    nxt = next_call(tool, last)
    if nxt is not None and nxt[0] == "cursor" and isinstance(answer, dict):
        return dict(answer, **{SELECTOR_CURSOR_KEY: nxt[1]})
    return answer
