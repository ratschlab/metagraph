"""The run facts a chapter's prose fills (round 4, after the re-verification of round 3, N-3): the chapter, family
and scenario values the prose names through placeholders and conditionals, computed from the record, so that a
later run's document never prints a figure of the run the prose was written for, and a sentence whose truth depends
on the run is printed only when the run makes it true.

Placeholders (str.format syntax, filled by `fill`): `{name}`; `{name:s}` a plural s when the number is not 1,
`{name:v}` the verb's s when it is 1, `{name:word}` a small number in words; `{x[R-pos-0027][dl][wall]}` a scenario
fact (the index syntax; `x` holds the chapter's or the family's scenarios, `dl`, `dr` and `mcp` the arms);
`{fam[R-budget][n_footnoted]}` a family's values from the chapter's prose. Conditionals: `[[flag ? text when true |
text when false]]`, `[[!flag ? ...]]`, `[[x[R-edge-0003][mcp][class] = refused_before_host ? ...]]`,
`[[manifest_version >= 5 ? ...]]` (a number at least the given one); nested ones are
resolved innermost first; the texts may hold placeholders; the false text is optional. A placeholder or flag the
values lack stays visible (`{name}`, `(UNRESOLVED CONDITION name)`) and is reported, never silently dropped.

The full document (the chapters T..D, 2026-10-15) adds: `manifest_version` at every level; per arm the first
execution's `time_limited` text and its `checks` by name ({outcome, detail, evidence}); the time limitations that
left an execution inconclusive, in three readings the chapters' prose uses (_time_limitations): the time cuts of a
walk (`n_time_cut_tests`, `time_cut_ids`, `time_cut_text`, `n_time_cut_untested_tests`, `time_cut_untested_ids`),
every such 200 (`n_time_limited`, `n_time_limited_tests`, `time_limited_ids`, `time_limited_text`) and the /pattern
ones (`n_time_limited_200`, `n_time_limited_200_tests`, `time_limited_200_ids`, `time_limited_200_text`);
`direct_all_pass`; and at the chapter level `unit` (every unit's unit_facts, read as `unit[<unit id>][...]`).

The v5 prose round (chapter R, 2026-10-15) adds, per family and per unit, the warm-ups' order against their targets
(_warm_ups: manifest v5's `warms`): `n_warm_ups`, `warm_up_targets_text`, `n_warm_up_arm_pairs`,
`n_warm_ups_after_target`, `warm_ups_all_before_target`, `warm_ups_after_target_text`, `n_warm_up_reruns` and
`warm_up_reruns_text`, read by R-warm's and R-budget's prose; after VERIFY-R5 the same order over the direct arms
alone (`n_warm_up_direct_pairs`, `n_warm_ups_after_target_direct`, `warm_ups_all_before_target_direct`,
`warm_ups_after_target_direct_text`) for R-budget, whose MCP calls the service refuses at submit;
`n_warm_up_reruns_without_target`; and per family `n_footnoted_late503` and `footnoted_all_late503`, so that
R-budget's prose claims 503 deadlines only for footnoted members that answered one. After VERIFY-R5-2 the 503 facts
add `late503_ids_text` ('a and b') and `late503_range_text` (one value when the earliest and the latest lateness
print alike), the warm-ups' order `warm_ups_after_target_direct_arms` (the direct arm, or 'both direct arms'), and per
family the stops that completed without a stop and without a 503 (`n_stops_complete`, `stops_complete_text`), so
that the coverage gap names the tests that answered a 503 and those that completed. After VERIFY-R5-3, per family
`n_stops` and `stops_all_complete` (every stop completed within its budget without a stop), so that the coverage gap says 'the stops completed'
only then and otherwise names those that did; and `n_late503_findings` (the 503 deadlines that are findings, beyond
max(1 s, the budget)) with `late503_findings_only_in_run` (the run's findings are exactly those), so that 'of these
503s' counts no other finding, such as a 200 past its budget. The run-wide placeholders
(sections.Context.placeholders) read `resolve_stops`: the stated stops of /resolve over the whole run on every arm
(`n_resolve_stops`, `resolve_stop_observed`, `resolve_stop_ids`), for the glossary's entry on the stated stop.

Every value is plain text (the prose is escaped afterwards; the inline markup `code` and *emphasis* is allowed),
a number or a bool; the ids are scenario ids as recorded.
"""

import os
import re
import string
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
if os.path.dirname(HERE) not in sys.path:
    sys.path.insert(0, os.path.dirname(HERE))

import blocks  # noqa: E402
import cells  # noqa: E402
import suite_common as sc  # noqa: E402
import units as units_mod  # noqa: E402
from cells import ARMS, ARM_SHORT, fmt_t  # noqa: E402

ARM_KEY = {"direct-local": "dl", "direct-remote": "dr", "mcp": "mcp"}
STANDING_WORD = {"pass": "pass", "fail": "fail", "error": "error", "finding": "finding", "inconclusive": "inconclusive",
                 "not_implemented_mcp": "not exercised", "service_by_design": "by design", "other": "other",
                 None: "no execution"}
HOST_WORD = {"cold": "host-cold", "warm": "host-warm", "aged": "host-aged",
             cells.WARM_AFTER_INTERRUPTED: "warm after an interrupted execution", "not_reached": "the host not reached",
             "unknown": "host state not recorded"}
WORDS = ["zero", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten", "eleven", "twelve"]
_COND = re.compile(r"\[\[((?:(?!\[\[|\]\]).)*)\]\]", re.S)
_NONE = ("", "0", "none", "–", "no", "False")
_MARKERS_RE = re.compile("[\ue000-\ue003]")


# -- filling ------------------------------------------------------------------------------------------------------

class _Blank(dict):
    """The facts of a test the run does not have (a cross reference, `xr`): every field reads as empty."""

    def __missing__(self, key):
        return _Blank()

    def __str__(self):
        return ""


class _Missing(object):
    def __init__(self, text):
        self.text = text

    def __str__(self):
        return self.text


class _Formatter(string.Formatter):
    """str.format over the values with the prose's extra specs; a missing field stays visible and is recorded."""

    def __init__(self, values):
        self.values = values
        self.unresolved = []

    def get_value(self, key, args, kwargs):
        if isinstance(key, str) and key in self.values:
            return self.values[key]
        raise KeyError(key)

    def get_field(self, field_name, args, kwargs):
        try:
            obj, first = string.Formatter.get_field(self, field_name, args, kwargs)
        except (KeyError, IndexError, TypeError, AttributeError):
            self.unresolved.append(field_name)
            return _Missing("{%s}" % field_name), field_name
        if obj is None:
            self.unresolved.append(field_name)
            return _Missing("{%s}" % field_name), field_name
        return obj, first

    def format_field(self, value, spec):
        if isinstance(value, _Missing):
            return str(value)
        n = cells._num(value)
        if n is None and spec in ("s", "v") and isinstance(value, str) and re.fullmatch(r"\d[\d,]*", value.strip()):
            # a count the closing values carry as text ('1'): its plural follows the number, not the type
            n = int(value.strip().replace(",", ""))
        if spec == "s":
            return "" if n == 1 else "s"
        if spec == "v":
            return "s" if n == 1 else ""
        if spec == "word":
            return number_word(value)
        if spec == "Word":
            w = number_word(value)
            return w[:1].upper() + w[1:]
        if spec == "":
            if isinstance(value, bool):
                return "yes" if value else "no"
            if isinstance(value, float):
                return ("%d" % value) if value == int(value) else ("%.1f" % value)
            return str(value)
        return string.Formatter.format_field(self, value, spec)


def number_word(value):
    """A small count in words (one ... twelve), larger ones as digits."""
    n = cells._num(value)
    if n is None:
        return str(value)
    if n == int(n) and 0 <= int(n) < len(WORDS):
        return WORDS[int(n)]
    return str(int(n)) if n == int(n) else str(n)


def _truthy(obj):
    if isinstance(obj, _Missing):
        return None
    if isinstance(obj, bool):
        return obj
    if isinstance(obj, (int, float)):
        return obj != 0
    if isinstance(obj, str):
        return obj.strip() not in _NONE
    return bool(obj)


def _splice(text, start, end, chosen):
    """`text` with [start, end) replaced by the chosen branch, the spaces around it as prose wants them: a branch
    that starts with a word after a word (or after a placeholder's closing brace) gets a space, a punctuation
    branch none; an empty branch leaves no double space and no space before a punctuation mark."""
    before, after = text[:start], text[end:]
    if chosen:
        if before and (before[-1].isalnum() or before[-1] in ")]}'’`") and chosen[0] not in ",.;:)!?":
            chosen = " " + chosen
        if after and (after[0].isalnum() or after[0] == "(" or (after[0] == "[" and not after.startswith("[["))):
            chosen = chosen + " "
        return before + chosen + after
    if before.endswith(" ") and after[:1] in (" ", ",", ".", ";", ":", ")"):
        before = before[:-1]
    return before + after


def fill(text, values, unresolved=None):
    """`text` with its conditionals resolved and its placeholders filled from `values`; the fields and flags the
    values lack are appended to `unresolved` (a list) and left visible."""
    fmt = _Formatter(values)
    text = str(text)
    # the conditionals, innermost first; a bracketed text without a `?` (a run list such as `[[0, 1]]`) is prose
    # and stays (so a literal `[[...]]` cannot sit inside a conditional's branch)
    pos = 0
    while True:
        m = _COND.search(text, pos)
        if not m:
            break
        body = m.group(1)
        head, sep, rest = body.partition("?")
        if not sep:
            pos = m.end()
            continue
        pos = 0
        yes, _, no = rest.partition("|")
        flag = head.strip()
        negate = flag.startswith("!")
        if negate:
            flag = flag[1:].strip()
        wanted = None
        at_least = ">=" in flag
        if at_least:
            # `[[manifest_version >= 5 ? ...]]`: a numeric comparison, so that a later manifest takes the branch of
            # the version that introduced the text (the v6 delta: v6 reads v5's sentences)
            flag, wanted = (p.strip() for p in flag.split(">=", 1))
        elif "=" in flag:
            flag, wanted = (p.strip() for p in flag.split("=", 1))
        obj, _ = fmt.get_field(flag, (), {})
        if isinstance(obj, _Missing):
            fmt.unresolved.pop()  # recorded under the condition instead
            fmt.unresolved.append("[[%s]]" % head.strip())
            text = text[:m.start()] + "(UNRESOLVED CONDITION %s)" % flag + text[m.end():]
            continue
        if at_least:
            have = cells._num(obj)
            try:
                want = float(wanted)
            except ValueError:
                want = None
            value = have is not None and want is not None and have >= want
        elif wanted is not None:
            # the printed value without the tint and attention markers (a phrase may carry them)
            value = (_MARKERS_RE.sub("", _Formatter(values).format_field(obj, "")) == wanted)
        else:
            value = bool(_truthy(obj))
        if negate:
            value = not value
        text = _splice(text, m.start(), m.end(), (yes if value else no).strip())
    try:
        out = fmt.vformat(text, (), {})
    except (ValueError, IndexError):
        out = text
        fmt.unresolved.append("(format error) " + text[:60])
    if unresolved is not None:
        unresolved.extend(fmt.unresolved)
    return out


# -- scenario facts ----------------------------------------------------------------------------------------------

def fmt_budget(ms):
    v = cells._num(ms)
    return "{:,} ms".format(int(v)) if v is not None else "no budget"


def _err(line):
    err = line.get("error") if isinstance(line.get("error"), dict) else {}
    return err.get("code"), err.get("text")


def _late_ms(run, sid, line):
    """How far past its budget a 503 deadline was delivered (client wall minus budget), or None."""
    code, _ = _err(line)
    budget = run.budget_ms(sid)
    wall = cells._num(line.get("wall_ms"))
    if line.get("http") == 503 and code == "deadline" and budget and wall is not None:
        return wall - budget
    return None


_ARM_HEADS = ("left", "right", "both")


def _clause(text, limit=160):
    """A check detail's evidence clause as plain text (the explanation after the first ': ' dropped). A walk
    check's detail names its arms ('left: the record followed for 793 of 1000 bp, stopped by time_budget; right:
    ...'): the first arm's clause up to ';' is kept with the arm's name; a 'not testable: ...' detail keeps the
    stop it names before the host's parenthesis (VERIFY-T M6: '(walk_in_record: left)' said nothing)."""
    text = str(text or "")
    head, rest = blocks._split_clause(text)
    if head.strip() in _ARM_HEADS and rest:
        seg, _ = blocks._split_clause(rest.split(";", 1)[0])
        return blocks._balanced(_cut("%s: %s" % (head.strip(), seg.strip()), limit))
    if head.strip() == "not testable" and rest:
        seg = rest.split(" (", 1)[0] if " (" in rest else blocks._split_clause(rest)[0]
        return blocks._balanced(_cut("not testable: %s" % seg.strip(), limit))
    return blocks._balanced(_cut(head, limit))


def _cut(text, limit):
    """`text` cut at `limit`, the execution times it quotes (a check's or the runner's text: cells.tint_quoted_times)
    carrying their tint as a marker that latex.esc sets."""
    import latex
    return cells.tint_quoted_times(latex.cut_text(str(text or ""), limit))


def _overridden(line, is_mcp):
    """The MCP runner's host-state rule overrode the checks (blocks._overridden): a scenario-level mark
    (suite_common.SCENARIO_LEVEL_MARK) is no host state, the checks decided it."""
    return blocks._overridden(line, is_mcp)


HOST_PHRASE = {"cold": "host-cold", "aged": "host-aged", cells.WARM_AFTER_INTERRUPTED: "warm after an interrupted execution"}


def _wall_text(line, arm):
    """The wall of an execution with its host state where that is not host-warm, on the direct arms ('1.9 s
    (host-cold)'; VERIFY-PLXM M5: a cold first call beside a warm one invites a timing reading the run does not
    isolate); the MCP arm's phrases name the cold and the warm call themselves."""
    text = cells.fmt_exec_text(line.get("wall_ms"))
    if arm != "mcp":
        word = HOST_PHRASE.get(cells.host_cache_of(line))
        if word:
            text += " (%s)" % word
    return text


def exec_phrase(run, sid, line, arm, wall=None):
    """One execution in a few plain words: the outcome, its qualifier and the time with its host state where not
    host-warm (`wall` overrides the time's text, so that repeats alike apart from their walls share one phrase)."""
    outcome = line.get("outcome")
    http = line.get("http")
    code, text = _err(line)
    wall = wall or _wall_text(line, arm)
    budget = run.budget_ms(sid)
    late = _late_ms(run, sid, line)
    if outcome == "pass":
        if http is not None and http >= 400:
            return "the expected refusal (http %s) in %s" % (http, wall)
        return "a pass in %s" % wall
    if outcome == "finding":
        if late is not None:
            return "a 503 deadline %s past the %s budget (a finding)" % (fmt_t(late), fmt_budget(budget))
        fc = [c for c in (line.get("checks") or []) if isinstance(c, dict) and c.get("outcome") == "finding"]
        return "a finding (%s)" % (_clause(fc[0].get("detail")) if fc else "a stated bound broken")
    if outcome == "inconclusive":
        if _overridden(line, arm == "mcp"):
            return "the service's %s at %s (a host state: inconclusive)" % (line["inconclusive_reason"], wall)
        if arm == "mcp" and line.get("inconclusive_reason") == sc.SCENARIO_LEVEL_MARK:
            return ("inconclusive at %s: every property check of this execution is marked as not shown via MCP (a "
                    "scenario-level mark), so nothing of the claim was tested through the tools" % wall)
        if late is not None:
            return "a 503 deadline %s past the %s budget (allowed: inconclusive)" % (fmt_t(late), fmt_budget(budget))
        if http == 503 and code == "deadline":
            return "a 503 deadline at %s (inconclusive)" % wall
        if http is not None and http != 200:
            return "http %s%s in %s (allowed: inconclusive)" % (http, (" " + code) if code else "", wall)
        checks = [c for c in (line.get("checks") or []) if isinstance(c, dict) and c.get("outcome") == "inconclusive"]
        if checks:
            return "inconclusive (%s: %s)" % (checks[0].get("check"), _clause(checks[0].get("detail"), 120))
        hidden = not_exposed_text(run, line)
        if hidden:
            # the MCP answer lacks the field the property check reads (the external review of the v6 document,
            # V6-302: the phrase names the cause, never a bare word)
            return "inconclusive (%s)" % hidden
        return "inconclusive"
    if outcome == "fail":
        failed = [c for c in (line.get("checks") or []) if isinstance(c, dict) and c.get("outcome") == "fail"]
        if failed:
            return "FAIL (%s: %s)" % (failed[0].get("check"), _cut(failed[0].get("detail"), 120))
        return "FAIL"
    if outcome == "error":
        return "ERROR (no HTTP answer%s)" % ((": " + _cut(text, 100)) if text else "")
    if outcome == "not_implemented_mcp":
        return cells.attn_marked("not exercised via MCP")
    if outcome == "service_by_design":
        return "by design (http %s%s)" % (http, (" " + code) if code else "")
    return str(cells.WORD.get(outcome, outcome))


def not_exposed_text(run, line):
    """'capabilities_doc not exposed via MCP: field not in the answer: features' for the property checks of an
    execution its answer did not expose (outcome not_exposed), '' when none (V6-302)."""
    items = [c for c in (line.get("checks") or []) if isinstance(c, dict) and c.get("outcome") == "not_exposed"
             and run.is_property(c.get("check"))]
    if not items:
        return ""
    return "; ".join("%s not exposed via MCP: %s" % (c.get("check"), _cut(c.get("detail"), 100)) for c in items)


def unshown_reason(run, sid, line, arm, short=False):
    """Why one execution did not show the test's property, without its time (the detail block's unverified line,
    V6-302): the scenario-level mark, a host state, a 503 deadline (allowed or late), another refusal, a failure, the
    property check that could not decide (with its evidence clause; `short`: the check alone, the exception index's
    grouping key) or that the MCP answer did not expose; '' for an execution that showed it."""
    outcome = line.get("outcome")
    code, _ = _err(line)
    http = line.get("http")
    if outcome == "finding":
        if _late_ms(run, sid, line) is not None:
            return "a 503 deadline past the budget (a finding, no answer)"
        props = [c for c in (line.get("checks") or []) if isinstance(c, dict) and run.is_property(c.get("check"))]
        if props and not any(c.get("outcome") == "pass" for c in props):
            return "a finding without a passing property check"
        return ""
    if outcome == "fail":
        failed = [c for c in (line.get("checks") or []) if isinstance(c, dict) and c.get("outcome") == "fail"]
        if failed and short:
            return "FAIL (%s)" % failed[0].get("check")
        return "FAIL (%s: %s)" % (failed[0].get("check"), _cut(failed[0].get("detail"), 100)) if failed else "FAIL"
    if outcome == "error":
        return "ERROR (no HTTP answer to check)"
    if outcome != "inconclusive":
        return ""
    if arm == "mcp" and line.get("inconclusive_reason") == sc.SCENARIO_LEVEL_MARK:
        return "every property check marked as not shown via MCP (a scenario-level mark)"
    if _overridden(line, arm == "mcp"):
        return "the service's %s (a host state)" % line["inconclusive_reason"]
    if http == 503 and code == "deadline":
        return "an allowed 503 deadline (no answer to check)"
    if http is not None and http != 200:
        return "http %s%s (allowed: no answer to check)" % (http, (" " + code) if code else "")
    items = [c for c in (line.get("checks") or []) if isinstance(c, dict) and c.get("outcome") == "inconclusive"]
    props = [c for c in items if run.is_property(c.get("check"))] or items
    if props:
        if short:
            return "%s could not decide" % props[0].get("check")
        # every property check the execution left undecided (the verification of the review's fixes, L6)
        return "; ".join("%s inconclusive: %s" % (c.get("check"), _clause(c.get("detail"), 110 if len(props) == 1
                                                                            else 80)) for c in props)
    if short:
        hidden = [c for c in (line.get("checks") or []) if isinstance(c, dict) and c.get("outcome") == "not_exposed"
                  and run.is_property(c.get("check"))]
        if hidden:
            fields = [str(c.get("detail") or "").split("field not in the answer: ", 1)[-1] for c in hidden]
            return "%s not exposed via MCP (the answer lacks %s)" % (
                ", ".join(str(c.get("check")) for c in hidden), ", ".join(_cut(f, 40) for f in fields))
    return not_exposed_text(run, line) or "no check could decide"


def unshown_reasons_text(run, sid, states):
    """Per access path the reasons the test's property was not shown there (unshown_reason), repeats alike
    grouped: 'mcp r0, r1: capabilities_doc not exposed via MCP: field not in the answer: features' (V6-302)."""
    parts = []
    for arm in ARMS:
        st = states[arm]
        if not st.present:
            continue
        by_reason = []
        marked = set(run.check_marks(sid)) if hasattr(run, "check_marks") else set()
        for l in st.lines:
            why = unshown_reason(run, sid, l, arm)
            if not why and l.get("outcome") == "pass":
                # a passing repeat whose property check the MCP answer did not expose, unless the manifest's marks
                # set that check aside (L6: v5's D-P-scope-0004 r1, r2 pass on deterministic alone)
                hidden = [c for c in (l.get("checks") or []) if isinstance(c, dict) and c.get("outcome") ==
                          "not_exposed" and run.is_property(c.get("check")) and c.get("check") not in marked]
                if hidden:
                    why = "; ".join("%s not exposed via MCP: %s" % (c.get("check"), _cut(c.get("detail"), 100))
                                    for c in hidden)
            if not why:
                continue
            for r in by_reason:
                if r[0] == why:
                    r[1].append(l.get("repeat"))
                    break
            else:
                by_reason.append([why, [l.get("repeat")]])
        for why, reps in by_reason:
            many = len(st.lines) > 1
            parts.append("%s%s: %s" % (arm, (" " + ", ".join("r%s" % r for r in reps)) if many else "", why))
    return "; ".join(parts)


def unshown_index_groups(run, sids, state_of):
    """The exception index's 'property not shown' entry (V6-302): the tests grouped by their recorded reason and
    the access paths it struck, [((access paths, reason), [ids])], e.g. (('mcp', 'capabilities_doc not exposed via
    MCP (the answer lacks features)'), ['C-caps-0003'])."""
    groups = []
    for sid in sids:
        states = state_of(sid)
        by_reason = []
        for arm in ARMS:
            st = states[arm]
            if not st.present:
                continue
            for l in st.lines:
                why = unshown_reason(run, sid, l, arm, short=True)
                if not why:
                    continue
                for r in by_reason:
                    if r[0] == why:
                        if arm not in r[1]:
                            r[1].append(arm)
                        break
                else:
                    by_reason.append([why, [arm]])
        for why, arms in by_reason:
            key = (_join(arms), why)
            for g in groups:
                if g[0] == key:
                    g[1].append(sid)
                    break
            else:
                groups.append([key, [sid]])
    return [(key, ids) for key, ids in groups]


def _status_text(line):
    code, _ = _err(line)
    http = line.get("http")
    if http is None:
        return ""
    return "http %s%s" % (http, (" " + code) if code else "")


def arm_facts(ctx, sid, arm):
    """The facts of one arm of a scenario: {standing, pass, present, wall, walls, wall_cold, wall_warm, server,
    backend, http, code (the first answered execution's status and error code; on an MCP arm the service decided,
    the first answered of the case's executions, the same the class is read from: blocks.class_lines), host,
    host_word (the first execution's host-cache state in words: host-warm, host-cold,
    warm after an interrupted execution), late, phrase, phrase_cold, phrase_warm, detail (the first property check's
    detail of the first execution), evidence (its clause before the explanation), class, class_word, clean,
    time_limited (the first execution's time limitation as the runner recorded it, '' when none), checks (the first
    execution's checks by name: {outcome, detail, evidence})}."""
    run = ctx.run
    st = ctx.state(sid)[arm]
    f = {"standing": STANDING_WORD.get(st.glyph_key, str(st.glyph_key)) if st.present else "no execution",
         "pass": bool(st.present and st.glyph_key == "pass"), "present": bool(st.present),
         "wall": "–", "walls": "–", "wall_cold": "–", "wall_warm": "–", "server": "–", "backend": "–",
         "http": "", "code": "", "host": "", "host_word": "", "late": "", "phrase": "", "phrase_cold": "",
         "phrase_warm": "", "detail": "", "evidence": "", "class": "", "class_word": "", "clean": bool(st.clean),
         "time_limited": "", "checks": {}, "n_outside_case": 0, "outside_case_phrase": "", "seed_refusal": ""}
    if not st.present:
        f["phrase"] = "no execution recorded although the run planned one" if st.missing_planned else "no execution"
        return f
    lines = st.lines
    first = lines[0]
    f["wall"] = cells.fmt_exec_text(first.get("wall_ms"))
    f["walls"] = " / ".join(cells.fmt_exec_text(l.get("wall_ms")) for l in lines)
    for l in lines:
        hc = cells.host_cache_of(l)
        if hc == "cold" and f["wall_cold"] == "–":
            f["wall_cold"] = cells.fmt_exec_text(l.get("wall_ms"))
        if hc == "warm" and f["wall_warm"] == "–":
            f["wall_warm"] = cells.fmt_exec_text(l.get("wall_ms"))
            srv = (l.get("server") or {}).get("elapsed_ms") if isinstance(l.get("server"), dict) else None
            f["server"] = cells.fmt_exec_text(srv)
            be = (l.get("mcp_phases") or {}).get("backend_ms") if isinstance(l.get("mcp_phases"), dict) else None
            f["backend"] = cells.fmt_exec_text(be)
    # on an MCP arm the status and the code come from the executions the class is read from (blocks.class_lines:
    # the case's executions on an arm the service decided), so that [mcp][class_word], [mcp][http] and [mcp][code]
    # describe one execution where the case covers one repeat only (VERIFY-FINAL R1, VERIFY-FIX2 L3)
    answered = next((l for l in (blocks.class_lines(st) if arm == "mcp" else lines) if l.get("http") is not None), None)
    if answered:
        f["http"] = str(answered.get("http"))
        f["code"] = _err(answered)[0] or ""
        # which seed refusal of SPEC-labeled-traversal-core section 6.1 the answer gave, read from its text (the
        # external review of the v6 document, V6-103): graph (step 2, a seed k-mer in no record of the graph: a plain
        # absence) or label (step 3, no seed label supports every k-mer); '' otherwise
        text = str(_err(answered)[1] or "")
        f["seed_refusal"] = ("graph" if "not fully present in the graph" in text else
                             "label" if "No seed label supports every k-mer" in text else "")
    f["host"] = str(cells.host_cache_of(first) or "").replace("_", " ")
    f["host_word"] = HOST_WORD.get(cells.host_cache_of(first), f["host"])
    lates = [_late_ms(run, sid, l) for l in lines]
    lates = [x for x in lates if x is not None]
    if lates:
        f["late"] = fmt_t(max(lates))
    f["time_limited"] = str(first.get("time_limited") or "")
    for c in first.get("checks") or []:
        if isinstance(c, dict) and c.get("check") and c["check"] not in f["checks"]:
            f["checks"][c["check"]] = {"outcome": str(c.get("outcome") or ""),
                                       "detail": blocks._balanced(_cut(c.get("detail"), 400)),
                                       "evidence": _clause(c.get("detail"), 200)}
    props = [c for c in (first.get("checks") or []) if isinstance(c, dict)
             and (c.get("property") is True or (c.get("property") is None and run.is_property(c.get("check"))))]
    if props:
        f["detail"] = blocks._balanced(_cut(props[0].get("detail"), 400))
        f["evidence"] = _clause(props[0].get("detail"), 200)
    if arm == "mcp":
        # the class of a by-design arm is read from the case's executions (blocks.class_lines; VERIFY-FINAL R1)
        cls, _ = blocks.mcp_class(blocks.class_lines(st), run.scenario(sid), st.agreement, run.cases)
        f["class"] = cls
        f["class_word"] = blocks.MCP_CLASS_WORD.get(cls, cls)
    if st.glyph_key == "not_implemented_mcp":
        f["phrase"] = cells.attn_marked("not exercised via MCP")
        if f["class"] == "refused_before_host" and answered:
            f["phrase"] += " (refused at submit: %s)" % _status_text(answered)
        return f
    if st.glyph_key == "service_by_design":
        covered = next((l for l in lines if l.get("http") is not None and cells.by_design_line(st, l)), None) or answered
        f["phrase"] = "by design (%s)" % (_status_text(covered) if covered else "no status recorded")
        if "u" in st.glyph_marks:
            f["phrase"] += ", the pair unverified (%s)" % blocks.unverified_reason(st.agreement, st)
        # the other repeat the case did not match: a separate host answer with its recorded outcome
        others = blocks.other_repeats(st)
        f["n_outside_case"] = len(others)
        f["outside_case_phrase"] = "; ".join("r%s outside the case %s" % (r, exec_phrase(run, sid, l, arm)) for r, l in others)
        if others:
            f["phrase"] += "; " + f["outside_case_phrase"]
        return f
    phrases = [exec_phrase(run, sid, l, arm) for l in lines]
    if arm == "mcp":
        cold = [p for l, p in zip(lines, phrases) if cells.host_cache_of(l) == "cold"]
        warm = [p for l, p in zip(lines, phrases) if cells.host_cache_of(l) == "warm"]
        f["phrase_cold"] = cold[0] if cold else ""
        f["phrase_warm"] = warm[0] if warm else ""
        if len(lines) == 2 and cold and warm:
            text = "the cold call %s, the warm call %s" % (cold[0], warm[0])
        elif len(lines) == 1:
            text = phrases[0]
        else:
            text = "; ".join("r%s %s" % (l.get("repeat"), p) for l, p in zip(lines, phrases))
    elif len(phrases) == 1:
        text = phrases[0]
    elif len(set(phrases)) == 1:
        text = ("both repeats " if len(phrases) == 2 else "every repeat ") + phrases[0]
    elif len(set(exec_phrase(run, sid, l, arm, wall="\x00") for l in lines)) == 1:
        # the repeats alike apart from their walls: one phrase, the walls in repeat order
        text = exec_phrase(run, sid, lines[0], arm, wall="\x00").replace(
            "\x00", " / ".join(_wall_text(l, arm) for l in lines)) + " (r%s)" % " / r".join(
            str(l.get("repeat")) for l in lines)
    else:
        text = "; ".join("r%s %s" % (l.get("repeat"), p) for l, p in zip(lines, phrases))
    if st.ne:
        text += " (DISAGREES with the reference access path)"
    elif "t" in st.glyph_marks:
        text += " (the two direct access paths differ by timing alone)"
    elif "s" in st.glyph_marks:
        text += " (another answer, the same check outcomes)"
    elif "n" in st.glyph_marks:
        text += " (not comparable with the reference access path)"
    f["phrase"] = text
    return f


def scenario_facts(ctx, sid):
    """The facts of one scenario: id, variant, title, unit, family, footnoted (not clean), the arms under dl, dr
    and mcp (arm_facts), and `phrase`, the arms' phrases joined."""
    meta = ctx.units.by_scenario.get(sid) or {}
    states = ctx.state(sid)
    f = {"id": sid, "variant": units_mod.variant_text(meta.get("variant") or ""),
         "title": str(ctx.run.scenario(sid).get("title") or ""),
         "unit": str(meta.get("unit_id") or ctx.units.unit_of(sid) or ""), "family": str(meta.get("family") or ""),
         "footnoted": not cells.is_clean(states)}
    for arm in ARMS:
        f[ARM_KEY[arm]] = arm_facts(ctx, sid, arm)
    f["phrase"] = "; ".join("%s %s" % (arm, f[ARM_KEY[arm]]["phrase"]) for arm in ARMS if f[ARM_KEY[arm]]["present"])
    # whether the deploy gate's policy selects the scenario, and on which access paths (READ-318: a chapter note
    # states the gate relation from the policy, never as a condition left to the reader)
    pol = ctx.run.policy if isinstance(ctx.run.policy, dict) else {}
    sel = pol.get("selection") if isinstance(pol.get("selection"), dict) else {}
    f["gate_selected"] = sid in set(str(x) for x in (sel.get("scenarios") or []))
    f["gate_arms"] = _join(str(a) for a in (pol.get("arms") or [])) or "no access path"
    return f


def _label(f):
    """'R-pos-0027 (2,000 bp, pfa)': the id with its variant; a variant's own parenthesis (the record key) joins the
    label's with a comma, never nested."""
    v = f["variant"]
    if not v:
        return f["id"]
    if v.endswith(")") and " (" in v:
        v = v[:-1].replace(" (", ", ", 1)
    return "%s (%s)" % (f["id"], v)


def _members_phrase(f, exercised_only=True):
    """A scenario's arms in words: 'sid (variant): direct-local ...; direct-remote ...; mcp ...' (an MCP arm standing
    not exercised or by design left out when `exercised_only`)."""
    parts = []
    for arm in ARMS:
        a = f[ARM_KEY[arm]]
        if not a["present"]:
            continue
        if exercised_only and arm == "mcp" and a["standing"] in ("not exercised", "by design") and not a["n_outside_case"]:
            continue
        parts.append("%s %s" % (arm, a["phrase"]))
    return "%s: %s" % (_label(f), "; ".join(parts) if parts else "no exercised access path")


def _join(items):
    """'a, b and c'."""
    items = list(items)
    if len(items) <= 1:
        return "".join(items)
    return ", ".join(items[:-1]) + " and " + items[-1]


def _ids_short(ids):
    """Scenario ids with a shared prefix once: 'REL-prefix_equals-0001, -0002 and -0003'."""
    ids = list(ids)
    if len(ids) <= 1:
        return _join(ids)
    prefix = ids[0].rsplit("-", 1)[0] + "-"
    if all(i.startswith(prefix) for i in ids):
        return _join([ids[0]] + ["-" + i[len(prefix):] for i in ids[1:]])
    return _join(ids)


def ids_ranges(ids):
    """Test ids with consecutive numbers as ranges, a shared prefix once per run of ids (the external review of the
    v6 document, V6-003): 'P-exact-0001 to -0008, P-scope-0001 to -0003, -0005, -0007 to -0012, P-edge-0004'; every
    id is still named, by its range."""
    groups = []
    for sid in ids:
        prefix, _, num = str(sid).rpartition("-")
        if not prefix or not num.isdigit():
            groups.append([str(sid), None, []])
            continue
        if groups and groups[-1][0] == prefix and groups[-1][1] is not None:
            groups[-1][2].append(num)
        else:
            groups.append([prefix, len(num), [num]])
    parts = []
    for prefix, width, nums in groups:
        if width is None:
            parts.append(prefix)
            continue
        runs = []
        for n in nums:
            if runs and int(n) == int(runs[-1][-1]) + 1:
                runs[-1].append(n)
            else:
                runs.append([n])
        texts = []
        for i, r in enumerate(runs):
            head = ("%s-%s" % (prefix, r[0])) if i == 0 else "-%s" % r[0]
            texts.append(head if len(r) == 1 else "%s to -%s" % (head, r[-1]) if len(r) > 2 else
                         "%s, -%s" % (head, r[-1]))
        parts.append(", ".join(texts))
    return ", ".join(parts)


# -- family and chapter facts -------------------------------------------------------------------------------------

def _standing_counts(ctx, sids, prefix=""):
    """n_pass_dl ... per arm and glyph over the scenarios' standing (by glyph, as the tables count)."""
    v = {}
    for arm in ARMS:
        k = ARM_KEY[arm]
        counts = {}
        for sid in sids:
            st = ctx.state(sid)[arm]
            key = st.glyph_key if st.present else None
            counts[key] = counts.get(key, 0) + 1
        v["%sn_pass_%s" % (prefix, k)] = counts.get("pass", 0)
        v["%sn_fail_%s" % (prefix, k)] = counts.get("fail", 0) + counts.get("error", 0) + counts.get("other", 0)
        v["%sn_finding_%s" % (prefix, k)] = counts.get("finding", 0)
        v["%sn_inconclusive_%s" % (prefix, k)] = counts.get("inconclusive", 0)
        v["%sn_missing_%s" % (prefix, k)] = counts.get(None, 0)
    mcp = [ctx.state(sid)["mcp"] for sid in sids]
    v[prefix + "n_mcp_not_exercised"] = sum(1 for st in mcp if st.present and st.glyph_key == "not_implemented_mcp")
    # the tests whose MCP access path does not stand ○ but one of whose executions had a property check set aside by
    # the manifest's property-level mark not implemented via MCP (ArmState.ni_marked): they count under not
    # implemented via MCP in the attention tags (attention_facts)
    prop = [sid for sid in sids if ctx.state(sid)["mcp"].ni_marked]
    v[prefix + "n_mcp_property_ni"] = len(prop)
    v[prefix + "mcp_property_ni_ids"] = _ids_short(prop) or "none"
    v[prefix + "n_mcp_by_design"] = sum(1 for st in mcp if st.present and st.glyph_key == "service_by_design")
    v[prefix + "n_mcp_by_design_unverified"] = sum(1 for st in mcp if st.present and st.glyph_key == "service_by_design"
                                                   and "u" in st.glyph_marks)
    v[prefix + "n_mcp_by_design_verified"] = v[prefix + "n_mcp_by_design"] - v[prefix + "n_mcp_by_design_unverified"]
    v[prefix + "mcp_by_design_all_verified"] = v[prefix + "n_mcp_by_design_unverified"] == 0
    # why a pair is unverified (blocks.unverified_cause): the direct side did not match the case, or the other MCP
    # repeat outside the case has no direct answer to compare with (VERIFY-QCND H3); a pair may count under both
    causes, texts = {}, []
    for sid in sids:
        st = ctx.state(sid)["mcp"]
        if st.present and st.glyph_key == "service_by_design" and "u" in st.glyph_marks:
            for c in (blocks.unverified_cause(st.agreement) or "direct").split("+"):
                causes.setdefault(c, []).append(sid)
            texts.append("%s (%s)" % (sid, blocks.unverified_reason(st.agreement, st)))
    for c in ("direct", "other_repeat"):
        v["%sn_mcp_by_design_unverified_%s" % (prefix, c)] = len(causes.get(c) or [])
        v["%smcp_by_design_unverified_%s_ids" % (prefix, c)] = ", ".join(causes.get(c) or []) or "none"
    v[prefix + "mcp_by_design_unverified_text"] = "; ".join(texts)
    return v


def _findings(ctx, sids):
    """(executions, tests, text) of the findings over the scenarios: 'sid 6.1 s late on direct-local and 7.5 s late
    on direct-remote; sid2 ...'."""
    run = ctx.run
    n_exec = 0
    per = []
    for sid in sids:
        items = []
        for arm in ARMS:
            for l in run.arm_lines(sid, arm):
                if l.get("outcome") != "finding":
                    continue
                n_exec += 1
                late = _late_ms(run, sid, l)
                if late is not None:
                    items.append("%s late on %s" % (fmt_t(late), arm))
                else:
                    fc = [c for c in (l.get("checks") or []) if isinstance(c, dict) and c.get("outcome") == "finding"]
                    items.append("%s on %s" % (_clause(fc[0].get("detail"), 100) if fc else "a stated bound broken", arm))
        if items:
            per.append("%s %s" % (sid, _join(items)))
    return n_exec, len(per), "; ".join(per)


def _late503(ctx, sids):
    """The 503 deadlines past a budget over the scenarios: {n_late503 (the executions), n_late503_tests, late503_ids
    ('a, b'), late503_ids_text ('a, b and c'), late503_min, late503_max, late503_range_text ('1.2 s to 25 s', or the one
    value when both print alike), late503_text, n_late503_findings (the executions among them that are findings: more
    than max(1 s, the budget) late, checks.budget_respected), and per arm n_late503_tests_<arm> and late503_<arm>_ids};
    `late503_id_list` is the ids as a list for the Python readers (family_facts), not for the prose."""
    run = ctx.run
    lates = []
    ids = []
    per = []
    n_findings = 0
    by_arm = dict((arm, []) for arm in ARMS)
    for sid in sids:
        items = []
        for arm in ARMS:
            for l in run.arm_lines(sid, arm):
                late = _late_ms(run, sid, l)
                if late is None:
                    continue
                lates.append(late)
                n_findings += 1 if l.get("outcome") == "finding" else 0
                items.append("%s on %s%s" % (fmt_t(late), arm, " (a finding)" if l.get("outcome") == "finding" else ""))
                if sid not in by_arm[arm]:
                    by_arm[arm].append(sid)
        if items:
            ids.append(sid)
            per.append("%s %s" % (sid, _join(items)))
    lo, hi = (fmt_t(min(lates)), fmt_t(max(lates))) if lates else ("–", "–")
    out = {"n_late503": len(lates), "n_late503_tests": len(ids), "late503_ids": ", ".join(ids),
           "late503_ids_text": _join(ids), "late503_id_list": list(ids), "late503_min": lo, "late503_max": hi,
           "late503_range_text": lo if lo == hi else "%s to %s" % (lo, hi), "late503_text": "; ".join(per),
           "n_late503_findings": n_findings}
    # the tests with a 503 deadline per arm (VERIFY-PLXM M4: the counts per direct arm)
    for arm in ARMS:
        out["n_late503_tests_%s" % ARM_KEY[arm]] = len(by_arm[arm])
        out["late503_%s_ids" % ARM_KEY[arm]] = ", ".join(by_arm[arm]) or "none"
    return out


def _stops(ctx, sids):
    """The stated stops answered (a 200 whose answer the checks read as stopped by time): {n_stopped,
    stop_observed, stop_observed_text}."""
    run = ctx.run
    per = []
    for sid in sids:
        for arm in ARMS:
            for l in run.arm_lines(sid, arm):
                tl = str(l.get("time_limited") or "")
                if l.get("http") == 200 and "stop" in tl:
                    per.append("%s on %s (%s, %s)" % (sid, arm, tl, cells.fmt_exec_text(l.get("wall_ms"))))
    return {"n_stopped": len(per), "stop_observed": bool(per), "stop_observed_text": "; ".join(per)}


def resolve_stops(run):
    """The stated stops of /resolve over the whole run, every scenario of the run's order whose route is /resolve on
    every arm (the glossary's 'stated stop, prefix', appendix D: a 200 whose recorded time limitation names a
    stop, checks.time_limited's 'resolve stopped by time'; a /traverse seed stopped by its time budget is not one):
    {n_resolve_stops (the executions), resolve_stop_observed, resolve_stop_ids (the scenarios in run order, 'none'
    when none)}. Read by the run-wide placeholders (sections.Context.placeholders), so that the glossary states
    from the record whether the run answered one, never as a fixed sentence."""
    n, ids = 0, []
    for sid in run.order:
        if "/resolve" not in str((run.scenario(sid) or {}).get("route") or ""):
            continue
        for arm in ARMS:
            for l in run.arm_lines(sid, arm):
                if l.get("http") == 200 and "stop" in str(l.get("time_limited") or ""):
                    n += 1
                    if sid not in ids:
                        ids.append(sid)
    return {"n_resolve_stops": n, "resolve_stop_observed": bool(n), "resolve_stop_ids": ", ".join(ids) or "none"}


_TIME_CUT_WORDS = ("time budget", "time_budget", "attempt")


def _arm_repeat(run, sid, arm, line):
    """'direct-local', or 'direct-local r1' where the arm recorded more than one execution of the scenario."""
    if len(run.arm_lines(sid, arm)) > 1:
        return "%s r%s" % (arm, line.get("repeat"))
    return arm


def _time_limitations(ctx, sids):
    """The 200s a time limitation left inconclusive (an execution answered 200, outcome inconclusive, whose
    `time_limited` is non-empty: a /pattern entry withheld or cut by its deadline, a seed stopped by its time budget
    or its attempt's bound, a resolve stopped by time); a 200 that passed although limited is not counted, since the
    claim it carries was shown. Three readings: every such execution ({n_time_limited, n_time_limited_tests,
    time_limited_ids, time_limited_text: 'sid on arm (the limitation)'}); the /pattern ones, whose limitation names a
    pattern ({n_time_limited_200, n_time_limited_200_tests, time_limited_200_ids, time_limited_200_text: 'sid wall on
    arm'}); and the time cuts of a walk, whose limitation names a time budget or the attempt's bound
    ({n_time_cut_tests, time_cut_ids, time_cut_text: 'sid on arm r<n>: the walk check's detail', and of those the
    tests no arm passed, n_time_cut_untested_tests and time_cut_untested_ids: the walk property rests on no evidence
    from the run})."""
    run = ctx.run
    all_ids, all_per, pat_ids, pat_per, pat_n, cut_ids, cut_per, untested = [], [], [], [], 0, [], [], []
    by_arm = dict((arm, {"all": [], "pat": []}) for arm in ARMS)
    for sid in sids:
        items, pat_items, cut_items = [], [], []
        for arm in ARMS:
            for l in run.arm_lines(sid, arm):
                tl = str(l.get("time_limited") or "")
                if not tl or l.get("http") != 200 or l.get("outcome") != "inconclusive":
                    continue
                where = _arm_repeat(run, sid, arm, l)
                items.append((where, tl))
                if sid not in by_arm[arm]["all"]:
                    by_arm[arm]["all"].append(sid)
                if "pattern" in tl:
                    pat_items.append("%s on %s" % (cells.fmt_exec_text(l.get("wall_ms")), where))
                    if sid not in by_arm[arm]["pat"]:
                        by_arm[arm]["pat"].append(sid)
                if any(w in tl for w in _TIME_CUT_WORDS):
                    checks = [c for c in (l.get("checks") or []) if isinstance(c, dict) and c.get("outcome") == "inconclusive"]
                    walk = next((c for c in checks if c.get("check") in ("walk_in_record", "walk_end")),
                                checks[0] if checks else None)
                    cut_items.append("on %s r%s: %s" % (arm, l.get("repeat"), blocks._balanced(_cut(
                        walk.get("detail"), 200)) if walk else tl))
        if items:
            all_ids.append(sid)
            by_text = {}
            for where, tl in items:
                by_text.setdefault(tl, []).append(where)
            all_per.append("; ".join("%s on %s (%s)" % (sid, _join(wheres), tl) for tl, wheres in by_text.items()))
        if pat_items:
            pat_ids.append(sid)
            pat_n += len(pat_items)
            pat_per.append("%s %s" % (sid, _join(pat_items)))
        if cut_items:
            cut_ids.append(sid)
            cut_per.append("%s %s" % (sid, _join(cut_items)))
            # untested: no execution on any arm and no repeat passed (VERIFY-T L6: per execution, not the arm's
            # worst over its repeats)
            if not any(l.get("outcome") == "pass" for arm in ARMS for l in run.arm_lines(sid, arm)):
                untested.append(sid)
    n_all = 0
    for sid in all_ids:
        for arm in ARMS:
            n_all += sum(1 for l in run.arm_lines(sid, arm) if str(l.get("time_limited") or "") and l.get("http") == 200
                         and l.get("outcome") == "inconclusive")
    out = {"n_time_limited": n_all, "n_time_limited_tests": len(all_ids), "time_limited_ids": ", ".join(all_ids),
           "time_limited_text": "; ".join(all_per),
           "n_time_limited_200": pat_n, "n_time_limited_200_tests": len(pat_ids),
           "time_limited_200_ids": ", ".join(pat_ids), "time_limited_200_text": "; ".join(pat_per),
           "n_time_cut_tests": len(cut_ids), "time_cut_ids": ", ".join(cut_ids), "time_cut_text": "; ".join(cut_per),
           "n_time_cut_untested_tests": len(untested), "time_cut_untested_ids": ", ".join(untested)}
    # the same per arm (VERIFY-PLXM M4: the withheld lists per direct arm, the MCP ones apart)
    for arm in ARMS:
        k = ARM_KEY[arm]
        out["n_time_limited_tests_%s" % k] = len(by_arm[arm]["all"])
        out["time_limited_%s_ids" % k] = ", ".join(by_arm[arm]["all"]) or "none"
        out["n_time_limited_200_tests_%s" % k] = len(by_arm[arm]["pat"])
        out["time_limited_200_%s_ids" % k] = ", ".join(by_arm[arm]["pat"]) or "none"
    # one text per arm: the tests a 503 deadline and the tests a time-limited 200 left inconclusive there
    # ('on direct-local six by a 503 deadline (…) and one by a 200 limited by time (…); on direct-remote …')
    parts = []
    for arm in ARMS:
        items = []
        d503 = [sid for sid in sids if any(_late_ms(run, sid, l) is not None for l in run.arm_lines(sid, arm))]
        if d503:
            items.append("%s by a 503 deadline (%s)" % (number_word(len(d503)), ", ".join(d503)))
        if by_arm[arm]["all"]:
            items.append("%s by a 200 limited by time (%s)" % (number_word(len(by_arm[arm]["all"])),
                                                              ", ".join(by_arm[arm]["all"])))
        if items:
            parts.append("%s %s" % ("via MCP" if arm == "mcp" else "on " + arm, " and ".join(items)))
    out["time_by_arm_text"] = "; ".join(parts)
    return out


_RESOURCE_STOP_WORDS = ("resource_limit", "max_steps", "max_live_paths", "max_paths", "max_output_bp")


def _budget_stops(ctx, sids):
    """_time_limitations and, with them, the walks a stated resource other than time stopped before they showed
    their property (an execution answered 200, outcome inconclusive, whose walk check reads 'stopped by
    resource_limit' or a path cap, the resource, its knob and its figures read from the answer's resource_stop
    through blocks.resource_stop_phrase; VERIFY-COLOURS-V5 item 2): {n_cap_stop_tests, cap_stop_ids,
    cap_stop_text: 'sid on arm r<n>: seed 0: the work budget (bounds.max_work_units, 100,000,000 units; used
    100,000,275) at 3,618 bp (left), ...; the walk check's clause'}, and every walk a budget stopped, time, memory
    or work, or the attempt's bound ({n_budget_stop_tests, budget_stop_ids, budget_stop_text}: the time cuts and the
    resource stops together, one entry per test)."""
    run = ctx.run
    out = _time_limitations(ctx, sids)
    cut_ids = [x for x in out["time_cut_ids"].split(", ") if x]
    cap_ids, cap_per = [], []
    for sid in sids:
        items = []
        for arm in ARMS:
            for l in run.arm_lines(sid, arm):
                if l.get("http") != 200 or l.get("outcome") != "inconclusive":
                    continue
                checks = [c for c in (l.get("checks") or []) if isinstance(c, dict) and c.get("outcome") == "inconclusive"
                          and c.get("check") in ("walk_in_record", "walk_end")]
                walk = next((c for c in checks if any(w in str(c.get("detail") or "") for w in _RESOURCE_STOP_WORDS)), None)
                if walk is None:
                    continue
                stop = blocks.resource_stop_phrase(run, l)
                clause = blocks._balanced(_cut(walk.get("detail"), 200))
                items.append("on %s r%s: %s" % (arm, l.get("repeat"), ("%s; %s" % (stop, clause)) if stop else clause))
        if items:
            cap_ids.append(sid)
            cap_per.append("%s %s" % (sid, _join(items)))
    budget_ids = list(cut_ids) + [sid for sid in cap_ids if sid not in cut_ids]
    by_sid = dict((t.split(" ", 1)[0], t) for t in out["time_cut_text"].split("; ") if t)
    for t in cap_per:
        sid = t.split(" ", 1)[0]
        by_sid[sid] = (by_sid[sid] + "; " + t.split(" ", 1)[1]) if sid in by_sid else t
    out.update({"n_cap_stop_tests": len(cap_ids), "cap_stop_ids": ", ".join(cap_ids), "cap_stop_text": "; ".join(cap_per),
                "n_budget_stop_tests": len(budget_ids), "budget_stop_ids": ", ".join(budget_ids),
                "budget_stop_text": "; ".join(by_sid[sid] for sid in budget_ids if sid in by_sid)})
    return out


def _relations(ctx, rids):
    """The relations' outcomes per arm: rel[rid][dl|dr|mcp], the counts per arm and outcome, and a text that
    groups the relations by their outcome pattern."""
    comp_rel = ctx.run.comp_relations()
    v = {"n_relations": len(rids), "rel": {}}
    for arm in ARMS:
        for o in ("pass", "fail", "inconclusive"):
            v["n_rel_%s_%s" % (o, ARM_KEY[arm])] = 0
    patterns = {}
    for rid in rids:
        arms_ = (comp_rel.get(rid) or {}).get("arms") or {}
        entry = {}
        for arm in ARMS:
            o = (arms_.get(arm) or {}).get("outcome") or "not exposed"
            entry[ARM_KEY[arm]] = o
            if o in ("pass", "fail", "inconclusive"):
                v["n_rel_%s_%s" % (o, ARM_KEY[arm])] += 1
        v["rel"][rid] = entry
        key = tuple(entry[ARM_KEY[arm]] for arm in ARMS)
        patterns.setdefault(key, []).append(rid)
    parts = []
    words = {"not_exposed": "not evaluated (not exposed)", "not exposed": "not evaluated (not exposed)",
             "not_evaluated": "not evaluated"}
    for key, ids in patterns.items():
        if len(set(key)) == 1:
            how = "%s on every access path" % words.get(key[0], key[0])
        else:
            by_outcome = {}
            for arm, o in zip(ARMS, key):
                by_outcome.setdefault(o, []).append(arm)
            how = ", ".join("%s on %s" % (words.get(o, o), _join(arms)) for o, arms in by_outcome.items())
        parts.append("%s %s" % (_ids_short(ids), how))
    v["relations_text"] = "; ".join(parts)
    return v


def _declared(ctx, sids, rids):
    """The owner's declared exemptions that name the scenarios or relations: counts and texts."""
    dec = ctx.run.declared if isinstance(ctx.run.declared, dict) else {}
    rel = [e for e in (dec.get("relations_excluded") or []) if isinstance(e, dict) and e.get("relation") in rids]
    exc = [e for e in (dec.get("exclusions") or []) if isinstance(e, dict) and e.get("scenario") in sids]
    by_arms = {}
    for e in rel:
        by_arms.setdefault(tuple(e.get("arms") or []), []).append(e["relation"])
    rel_text = "; ".join("%s on %s" % (_ids_short(ids), _join(arms) if arms else "every access path")
                         for arms, ids in by_arms.items())
    exc_text = "; ".join("%s on %s" % (e["scenario"], _join(e.get("arms") or []) or "every access path") for e in exc)
    return {"n_declared_relations": len(rel), "declared_relations_text": rel_text,
            "n_declared_exclusions": len(exc), "declared_exclusions_text": exc_text}


def _expected_runs(ctx, sids):
    """'variant: [0, 71) and [100, 170); ...' from the members' resolve_runs expectations that name runs."""
    per = []
    for sid in sids:
        runs = []
        for e in ctx.run.scenario(sid).get("expect") or []:
            if isinstance(e, dict) and e.get("check") == "resolve_runs" and (e.get("args") or {}).get("runs"):
                for r in e["args"]["runs"]:
                    if isinstance(r, (list, tuple)) and len(r) == 2:
                        runs.append("[%s, %s)" % (r[0], r[1]))
        if runs:
            variant = (ctx.units.by_scenario.get(sid) or {}).get("variant") or sid
            per.append("%s %s" % (variant, _join(runs)))
    return "; ".join(per)


def _gate_blocking(ctx, sids):
    gate = ctx.run.gate if isinstance(ctx.run.gate, dict) else {}
    ids = []
    for p in gate.get("problems") or []:
        if isinstance(p, dict) and p.get("kind") == "outcome" and p.get("scenario") in sids and p["scenario"] not in ids:
            ids.append(p["scenario"])
    return ids


def _cold_counts(ctx, sids):
    v = {}
    for arm in ("direct-local", "direct-remote"):
        v["n_cold_%s" % ARM_KEY[arm]] = sum(1 for sid in sids for l in ctx.run.arm_lines(sid, arm)
                                             if cells.host_cache_of(l) == "cold")
    return v


def _mcp_class_values(x, members):
    """The MCP classes over the members as the cells derive them (blocks.mcp_class): `n_mcp_<class>`,
    `mcp_<class>_ids`, `mcp_<class>_variants`, `mcp_<class>_codes` (the service's codes, as code spans);
    `mcp_classes_text` ('6 refused before host execution, 1 answered'); `mcp_all_refused_before_host` and
    `mcp_all_not_exercised` (every member), for a prose sentence that otherwise asserts the run."""
    v = {}
    by_class = {}
    for sid in members:
        by_class.setdefault(x[sid]["mcp"]["class"] or "none", []).append(sid)
    parts = []
    for cls in blocks.MCP_CLASS_WORD:
        ids = by_class.get(cls) or []
        v["n_mcp_%s" % cls] = len(ids)
        v["mcp_%s_ids" % cls] = ", ".join(ids) or "none"
        v["mcp_%s_variants" % cls] = ", ".join(x[s]["variant"] or s for s in ids) or "none"
        codes = sorted(set(x[s]["mcp"]["code"] for s in ids if x[s]["mcp"]["code"]))
        v["mcp_%s_codes" % cls] = ", ".join("`%s`" % c for c in codes) or "none"
        if ids:
            # a class whose members stand not exercised via MCP is a not-implemented count, red wherever printed
            # (the owner's decision of 2026-10-16); a class some of whose members do names those apart
            part = "%d %s" % (len(ids), blocks.MCP_CLASS_WORD.get(cls, cls))
            ni = sum(1 for sid in ids if x[sid]["mcp"]["standing"] == "not exercised")
            if ni == len(ids):
                part = cells.attn_marked(part)
            elif ni:
                part += " (%s)" % cells.attn_marked("%d of them not exercised" % ni)
            parts.append(part)
    v["mcp_classes_text"] = ", ".join(parts) or "no MCP execution"
    n = len(members)
    v["mcp_all_refused_before_host"] = bool(n) and v.get("n_mcp_refused_before_host", 0) == n
    v["mcp_all_not_exercised"] = bool(n) and all(x[s]["mcp"]["standing"] == "not exercised" for s in members)
    return v


def _started_s(line):
    """An execution's start in seconds since the epoch, or None (the merged lines carry started_utc)."""
    try:
        return sc.parse_utc(line.get("started_utc"))
    except Exception:
        return None


def _warm_ups(ctx, sids):
    """The warm-ups among the scenarios (manifest v5: `warms` names the scenario a warm-up precedes, which the build
    places right after it, build_manifest.place_warm_ups, so that the target's first execution meets rows the host
    has just read) against the run's own order: on each arm that ran both, whether the warm-up's first execution
    started before the target's first. {n_warm_ups; warm_up_targets_text ('R-warm-0001 (T-repeat-0001); ...');
    n_warm_up_arm_pairs (the (warm-up, arm) pairs compared); n_warm_ups_after_target (the pairs where the warm-up
    did not start before its target); warm_ups_all_before_target (every compared pair in order, at least one
    compared); warm_ups_after_target_text ('sid on arm (started 12 s after its target T-...)', a tie 'started with
    its target'); the same four over the direct arms alone (n_warm_up_direct_pairs, n_warm_ups_after_target_direct,
    warm_ups_all_before_target_direct, warm_ups_after_target_direct_text), for a target whose MCP call the service
    refuses at submit, where the order on that arm warms nothing (VERIFY-R5 L4), and warm_ups_after_target_direct_arms
    (the one direct arm those pairs are on, or 'both direct arms'; VERIFY-R5-2 issue 4); n_warm_up_reruns (the warm-ups
    among the scenarios run again for a resume, run.json reruns_for_resume: timing evidence, never an outcome);
    n_warm_up_reruns_without_target (those no execution of the target followed on that arm, target_execution
    null); warm_up_reruns_text ('sid on arm before T-repeat-0001', or 'sid on arm (no execution of its target T-...
    followed)')}."""
    run = ctx.run
    warm = [(sid, run.scenario(sid).get("warms")) for sid in sids if run.scenario(sid).get("warms")]
    pairs, after, per = 0, [], []
    direct_pairs, direct_after, direct_per, direct_arms = 0, [], [], []
    for sid, target in warm:
        for arm in ARMS:
            a = [t for t in (_started_s(l) for l in run.arm_lines(sid, arm)) if t is not None]
            b = [t for t in (_started_s(l) for l in run.arm_lines(target, arm)) if t is not None]
            if not a or not b:
                continue
            direct = arm != "mcp"
            pairs += 1
            direct_pairs += 1 if direct else 0
            if min(a) >= min(b):
                when = "started with" if min(a) == min(b) else "started %s after" % fmt_t((min(a) - min(b)) * 1000.0)
                item = "%s on %s (%s its target %s)" % (sid, arm, when, target)
                after.append(sid)
                per.append(item)
                if direct:
                    direct_after.append(sid)
                    direct_per.append(item)
                    if arm not in direct_arms:
                        direct_arms.append(arm)
    reruns = [r for r in run.reruns if r.get("scenario") in [sid for sid, _ in warm]]

    def rerun_text(r):
        if r.get("target_execution"):
            return "%s on %s before %s" % (r.get("scenario"), r.get("arm"), r.get("rerun_target"))
        return "%s on %s (no execution of its target %s followed)" % (r.get("scenario"), r.get("arm"), r.get("rerun_target"))

    return {"n_warm_ups": len(warm), "warm_up_targets_text": "; ".join("%s (%s)" % (sid, t) for sid, t in warm),
            "n_warm_up_arm_pairs": pairs, "n_warm_ups_after_target": len(after),
            "warm_ups_all_before_target": bool(pairs) and not after,
            "warm_ups_after_target_text": "; ".join(per),
            "n_warm_up_direct_pairs": direct_pairs, "n_warm_ups_after_target_direct": len(direct_after),
            "warm_ups_all_before_target_direct": bool(direct_pairs) and not direct_after,
            "warm_ups_after_target_direct_text": "; ".join(direct_per),
            "warm_ups_after_target_direct_arms": "both direct access paths" if len(direct_arms) > 1 else "".join(direct_arms),
            "n_warm_up_reruns": len(reruns),
            "n_warm_up_reruns_without_target": sum(1 for r in reruns if not r.get("target_execution")),
            "warm_up_reruns_text": "; ".join(rerun_text(r) for r in reruns)}


def _expects_prefix(scenario):
    """Whether the scenario is a stop: its resolve_runs expects `prefix_if_stopped`."""
    return any(isinstance(e, dict) and (e.get("args") or {}).get("prefix_if_stopped")
               for e in (scenario or {}).get("expect") or [])


def _completed_without_stop(run, sid):
    """Whether every direct execution of the scenario (at least one) answered a complete 200: a pass or an inconclusive
    outcome, no stop in its time limitation."""
    lines = [l for arm in ARMS if arm != "mcp" for l in run.arm_lines(sid, arm)]
    return bool(lines) and all(l.get("http") == 200 and l.get("outcome") in ("pass", "inconclusive")
                               and "stop" not in str(l.get("time_limited") or "") for l in lines)


def unit_facts(ctx, uid, x=None):
    """A unit's values for the line under its table (`mcp_notes.units`) and for the chapter prose's `unit[<id>]`:
    n_tests, the members' facts under `x` (the scenario facts already computed may be passed), the counts by
    standing, the MCP classes, the time limitations and the warm-ups' order (_warm_ups) over its members, as
    family_facts computes them for a family, `direct_all_pass` and `manifest_version`."""
    u = ctx.units.units.get(uid) or {}
    members = list(u.get("members") or [])
    x = dict((sid, (x or {}).get(sid) or scenario_facts(ctx, sid)) for sid in members)
    v = {"id": uid, "n_tests": len(members), "members": ", ".join(members), "x": x}
    v.update(_mcp_class_values(x, members))
    v.update(_standing_counts(ctx, members))
    v.update(_budget_stops(ctx, members))
    v.update(_warm_ups(ctx, members))
    v["n_direct_pass_both"] = _direct_pass_both(x, members)
    v["direct_all_pass"] = bool(members) and v["n_direct_pass_both"] == len(members)
    v["manifest_version"] = _manifest_version(ctx)
    return v


def _manifest_version(ctx):
    return (ctx.run.manifest or {}).get("manifest_version")


def _direct_pass_both(x, members):
    """The members that pass on direct-local and on direct-remote where an execution exists."""
    return sum(1 for sid in members if x[sid]["dl"]["pass"] and (x[sid]["dr"]["pass"] or not x[sid]["dr"]["present"]))


def family_facts(ctx, fid):
    """A family's values for its prose (the chapter's prose reads them under `fam[<fid>]`): the counts by standing,
    the MCP classes as the cells derive them (`n_mcp_<class>`, `mcp_<class>_ids`, `mcp_<class>_variants`,
    `mcp_<class>_codes`, `mcp_classes_text`, `mcp_all_refused_before_host`), the clean and footnoted members with
    their phrases (`clean_text`, `footnoted_text`), every member's phrase, the findings, the 503 deadlines past
    their budgets (and whether every footnoted member answered one: `n_footnoted_late503`,
    `footnoted_all_late503`; whether the run's findings are exactly the 503 deadlines that are findings:
    `late503_findings_only_in_run`), the stated stops, the stops (`n_stops`) and those that completed without a stop
    and without a 503 (`n_stops_complete`, `stops_complete_text`, `stops_all_complete`), the time limitations
    (_time_limitations), the relations, the owner's declared
    exemptions, the warm-ups' order against their targets (_warm_ups, manifest v5), the expected runs, the gate's
    outcome problems among the members, `direct_all_pass`, `manifest_version` and `x` (the members' scenario facts)."""
    f = ctx.units.families.get(fid) or {}
    members = [sid for uid in f.get("units") or [] for sid in ctx.units.units[uid]["members"]]
    rids = [rid for uid in f.get("units") or [] for rid in (ctx.units.units[uid].get("relations") or [])]
    prose = ((ctx.prose.get(f.get("chapter")) or {}).get("families") or {}).get(fid) or {}
    v = {"id": fid, "title": str(prose.get("title") or f.get("title") or fid), "n_tests": len(members),
         "n_units": len(f.get("units") or []), "members": ", ".join(members)}
    x = dict((sid, scenario_facts(ctx, sid)) for sid in members)
    v["x"] = x
    # the scenario facts of another family's test the prose names as `xr[<id>]` (a cross reference: R-rc reads
    # T-strand's observed refusal, the verification of the review's fixes, L1), kept apart from the members' `x`
    # a test another run lacks reads blank (every field empty), so that a condition on it takes its last branch
    v["xr"] = dict((sid, scenario_facts(ctx, sid) if sid in (ctx.run.scenarios or {}) else _Blank())
                   for sid in sorted(set(re.findall(r"xr\[([A-Za-z0-9_.-]+)\]", str(prose)))))
    v.update(_mcp_class_values(x, members))
    v.update(_standing_counts(ctx, members))
    footnoted = [sid for sid in members if x[sid]["footnoted"]]
    clean = [sid for sid in members if not x[sid]["footnoted"]]
    v["n_clean"] = len(clean)
    v["n_footnoted"] = len(footnoted)
    v["footnoted_ids"] = ", ".join(footnoted)
    v["footnoted_text"] = "; ".join(_members_phrase(x[sid]) for sid in footnoted)
    v["clean_ids"] = ", ".join(clean)
    v["clean_text"] = "; ".join(_members_phrase(x[sid]) for sid in clean)
    v["members_text"] = "; ".join(_members_phrase(x[sid]) for sid in members)
    v["n_direct_pass_both"] = _direct_pass_both(x, members)
    v["direct_all_pass"] = bool(members) and v["n_direct_pass_both"] == len(members)
    v["manifest_version"] = _manifest_version(ctx)
    v["direct_failed"] = any(x[sid][k]["standing"] in ("fail", "error", "other") for sid in members for k in ("dl", "dr")) \
        or any(ctx.state(sid)["direct-remote"].ne for sid in members)
    n_exec, n_tests, text = _findings(ctx, members)
    v["n_finding_executions"], v["n_finding_tests"], v["findings_text"] = n_exec, n_tests, text
    v.update(_late503(ctx, members))
    # whether every footnoted member answered a 503 deadline past its budget (VERIFY-R5 M1: R-budget's prose claims
    # the coverage gap's 503 sentence only then; a stop that completed within its budget is footnoted without one)
    late_ids = v["late503_id_list"]
    v["n_footnoted_late503"] = sum(1 for sid in footnoted if sid in late_ids)
    v["footnoted_all_late503"] = bool(footnoted) and v["n_footnoted_late503"] == len(footnoted)
    v.update(_stops(ctx, members))
    v.update(_budget_stops(ctx, members))
    v.update(_cold_counts(ctx, members))
    v.update(_relations(ctx, rids))
    v.update(_declared(ctx, members, rids))
    v.update(_warm_ups(ctx, members))
    v["expected_runs_text"] = _expected_runs(ctx, members)
    blocking = _gate_blocking(ctx, members)
    v["n_gate_blocking"] = len(blocking)
    v["gate_blocking_ids"] = ", ".join(blocking)
    stops = [sid for sid in members if _expects_prefix(ctx.run.scenario(sid))]
    v["stop_ids_text"] = ", ".join("%s under %s" % (sid, fmt_budget(ctx.run.budget_ms(sid))) for sid in stops
                                   if ctx.run.budget_ms(sid))
    # the stops that completed: every direct execution a complete 200 and no 503 deadline on any arm (VERIFY-R5-2
    # issue 1: the coverage gap says so of a stop even when another member answered a 503)
    complete = [sid for sid in stops if sid not in late_ids and _completed_without_stop(ctx.run, sid)]
    v["n_stops"] = len(stops)
    v["n_stops_complete"] = len(complete)
    v["stops_complete_text"] = _join(complete)
    v["stops_all_complete"] = bool(stops) and len(complete) == len(stops)
    # whether the run's findings are exactly this family's 503 deadlines beyond max(1 s, the budget) (VERIFY-R5-3
    # issue 1: the coverage gap's 'of these 503s ... the run's only findings' counts no other finding, such as a
    # 200 past its budget, and claims the run's findings only when they are these)
    totals = ctx.run.summary().get("totals") or {}
    run_findings = sum((totals.get(a) or {}).get("outcomes", {}).get("finding", 0) for a in ARMS)
    v["late503_findings_only_in_run"] = bool(v["n_late503_findings"]) and v["n_late503_findings"] == run_findings
    return v


def claim_items(fp):
    """A family's claims as [(text, members or None)]: the prose's `claim` (a string over `claim_members`, or a
    list of strings or of {claim, members} each over its own members, so that a family of two units states one
    claim per unit: VERIFY-PLXM H1, M-protein's residue claim apart from its cap claim)."""
    claim = fp.get("claim")
    if not claim:
        return []
    if isinstance(claim, list):
        out = []
        for item in claim:
            if isinstance(item, dict) and item.get("claim"):
                out.append((str(item["claim"]), item.get("members")))
            elif isinstance(item, str) and item:
                out.append((item, None))
        return out
    return [(str(claim), fp.get("claim_members"))]


def _claims(ctx, ch, fam):
    """The families' claims (claim_items: the prose's `claim`, over its members or every member) sorted by what
    the direct arms showed: held (every claim member passes on both direct arms), held in part (some), unsupported
    (none, or a failure). A claim's text is filled on the family's values (a conditional on `manifest_version`
    for a property only one manifest's member tests: VERIFY-T H4)."""
    prose = (ctx.prose.get(ch) or {}).get("families") or {}
    held, unsupported = [], []
    items = []  # the same claims as records, for the chapter's evidence map (READ-101): sections.claims_map
    for fid in ctx.units.chapter_families(ch):
        fp = prose.get(fid) or {}
        fv = fam[fid]
        for claim, wanted in claim_items(fp):
            _claim(ctx, fid, fv, fill(claim, fv, ctx.unresolved), wanted, held, unsupported, items)
    return held, unsupported, items


def _claim(ctx, fid, fv, claim, wanted, held, unsupported, items=None):
    """One claim sorted into `held` or `unsupported` by the direct access paths' standing on its members (`wanted`,
    else every member of the family); `items` collects the record {family, claim, n, m, status, words}: status
    held (every member passes on both direct paths), part (some), none (no member; `words` the standings seen) or
    failed (a failure, an error or a disagreement among the members)."""
    members = [sid for sid in (wanted or []) if sid in fv["x"]] or list(fv["x"])
    x = fv["x"]
    n = sum(1 for sid in members if x[sid]["dl"]["pass"] and (x[sid]["dr"]["pass"] or not x[sid]["dr"]["present"]))
    failed = any(x[sid][k]["standing"] in ("fail", "error", "other") for sid in members for k in ("dl", "dr")) \
        or any(ctx.state(sid)["direct-remote"].ne for sid in members)
    words = []
    for sid in members:
        for k in ("dl", "dr"):
            w = x[sid][k]["standing"]
            if x[sid][k]["present"] and w not in words:
                words.append(w)
    if failed:
        unsupported.append("%s (%s: %d of %d, with a failure or a disagreement; I.1)" % (claim, fid, n, len(members)))
        status = "failed"
    elif n == len(members):
        held.append("%s (%s)" % (claim, fid))
        status = "held"
    elif n:
        held.append("%s (%s, %d of %d)" % (claim, fid, n, len(members)))
        status = "part"
    else:
        unsupported.append("%s (%s: %s on the direct access paths)" % (claim, fid, " / ".join(words) or "no execution"))
        status = "none"
    if items is not None:
        items.append({"family": fid, "claim": claim, "n": n, "m": len(members), "status": status,
                      "words": " / ".join(words) or "no execution"})


MCP_FULL_DEF = "every property check evaluated on both"


def mcp_full_definition(ctx, marked):
    """What 'fully checked' covers in a chapter's MCP sentence (V6-201): every property check evaluated on both, or,
    where agreeing pairs omit a check the manifest's property-level marks set aside, that exception with their count
    by family; the marked check is not compared and stays unverified via MCP."""
    if not marked:
        return MCP_FULL_DEF
    run = getattr(ctx, "run", None)
    bd = [sid for sid in marked if hasattr(run, "check_marks")
          and any(m.get("mark") == "by_design" for m in run.check_marks(sid).values())]
    ni = [sid for sid in marked if sid not in bd]
    text = MCP_FULL_DEF
    if ni:
        text = _full_def_ni(text, ni)
    if bd:
        # a check set aside as a documented service difference, told apart from the not-implemented kind (RV4)
        named = ["%s's %s" % (sid, ", ".join(c for c, m in run.check_marks(sid).items() if m.get("mark") == "by_design"))
                 for sid in bd[:3]]
        text += "%s except %s, which appendix G lists as set aside service-side by design" % (
            ";" if ni else ",", _join(named) + (" and %d more" % (len(bd) - 3) if len(bd) > 3 else ""))
    return text


def _full_def_ni(base, marked):
    """'fully checked' with the exception of the not-implemented property-level marks (V6-201), by family when many."""
    if len(marked) <= 3:
        which = _ids_short(marked)
    else:
        fams = []
        for sid in marked:
            fid = sid.rsplit("-", 1)[0]
            for f in fams:
                if f[0] == fid:
                    f[1] += 1
                    break
            else:
                fams.append([fid, 1])
        which = ", ".join("%s %d" % (f, k) for f, k in fams)
    return ("%s, except a check a property-level mark of the manifest sets aside, not compared and unverified via MCP "
            "(on %s of them: %s)" % (base, number_word(len(marked)), which))


def chapter_facts(ctx, ch):
    """A chapter's values for its `summary` and `coverage_gap` prose: the counts by standing and by agreement, the
    timing of the summary's sources, the findings, the time limitations, the claims held and unsupported, the
    footnoted tests, the gate, `fam` (every family's values), `unit` (every unit's values) and `x` (every scenario's
    facts)."""
    run = ctx.run
    fams = ctx.units.chapter_families(ch)
    sids = ctx.units.chapter_scenarios(ch)
    fn = ctx.footnotes.get(ch)
    v = {"chapter": ch, "n_tests": len(sids), "n_families": len(fams),
         "n_units": sum(len(ctx.units.families[f]["units"]) for f in fams),
         "n_detail": len(fn.order) if fn else 0}
    fam = dict((fid, family_facts(ctx, fid)) for fid in fams)
    v["fam"] = fam
    x = {}
    for fid in fams:
        x.update(fam[fid]["x"])
    v["x"] = x
    v["unit"] = dict((uid, unit_facts(ctx, uid, x)) for fid in fams for uid in ctx.units.families[fid]["units"])
    v.update(_standing_counts(ctx, sids))
    v["n_direct_pass_both"] = _direct_pass_both(x, sids)
    v["direct_all_pass"] = bool(sids) and v["n_direct_pass_both"] == len(sids)
    v["manifest_version"] = _manifest_version(ctx)
    # the agreement classes, as the chapter's summary block counts them
    direct, mcp = {}, {}
    full, verified, unverified = 0, 0, 0
    fewer = []  # the fully checked pairs whose comparison left a repeat out (a host-state outcome)
    full_marked = []  # the fully checked pairs with a property check a property-level mark sets aside (V6-201)
    for sid in sids:
        e = run.comp_entry(sid)
        dc = (e.get("direct_agreement") or {}).get("class")
        mc = (e.get("mcp_agreement") or {}).get("class")
        if dc:
            direct[dc] = direct.get(dc, 0) + 1
        if mc:
            mcp[mc] = mcp.get(mc, 0) + 1
            if mc == "service_by_design":
                if (e.get("mcp_agreement") or {}).get("verified") is False:
                    unverified += 1
                else:
                    verified += 1
        d = (run.depth or {}).get(sid)
        if d and d.get("class") == "full":
            full += 1
            if d.get("repeats_excluded"):
                fewer.append(sid)
            if hasattr(run, "check_marks") and run.check_marks(sid):
                full_marked.append(sid)
    for key in ("identical", "same_checks", "timing", "not_comparable", "disagree"):
        v["n_direct_%s" % key] = direct.get(key, 0)
    words = {"identical": "identical pair{s}", "timing": "differing by a time limitation alone",
             "same_checks": "with the same check outcomes on different answers", "not_comparable": "not comparable",
             "disagree": "disagreeing"}
    v["direct_agreement_text"] = ", ".join(
        "%d %s" % (direct[k], words[k].replace("{s}", "" if direct[k] == 1 else "s"))
        for k in ("identical", "timing", "same_checks", "not_comparable", "disagree") if direct.get(k)) or "no pair"
    v["n_direct_failed_tests"] = sum(1 for sid in sids if any(x[sid][k]["standing"] in ("fail", "error", "other")
                                                                 for k in ("dl", "dr")) or ctx.state(sid)["direct-remote"].ne)
    v["direct_clean"] = v["n_direct_failed_tests"] == 0
    v["n_mcp_agree"] = mcp.get("agree", 0)
    v["n_mcp_full"] = full
    v["mcp_all_full"] = full == mcp.get("agree", 0)
    # fully checked holds over the property checks the manifest's property-level marks do not set aside (the
    # external review of the v6 document, V6-201): the pairs that omit a marked check, counted by family, so that
    # the chapter's definition never claims the marked check was evaluated via MCP
    v["n_mcp_full_marked"] = len(full_marked)
    v["mcp_full_def"] = mcp_full_definition(ctx, full_marked)
    # the fully checked pairs compared on fewer repeats than recorded, each with the repeats compared and the
    # excluded one's reason (reader.depth_over_compared_repeats; the owner's second review, point 2)
    v["n_mcp_full_fewer_repeats"] = len(fewer)
    v["mcp_full_fewer_repeats_text"] = "; ".join(
        "%s on %s, %s" % (sid, _join("r%s" % r for r in (run.depth[sid].get("repeats_compared") or [])),
                          _join("r%s %s" % (e.get("repeat"), e.get("reason"))
                                for e in run.depth[sid].get("repeats_excluded") or [] if isinstance(e, dict)))
        for sid in fewer)
    v["n_mcp_disagree"] = mcp.get("disagree", 0)
    v["n_mcp_not_comparable"] = mcp.get("not_comparable", 0)
    v["n_mcp_not_implemented"] = mcp.get("not_implemented_mcp", 0)
    v["n_mcp_by_design_class"] = mcp.get("service_by_design", 0)
    # the not-exercised tests by family, with the classes where a family mixes them
    parts = []
    for fid in fams:
        fv = fam[fid]
        n = fv["n_mcp_not_exercised"]
        if not n:
            continue
        classes = {}
        for sid in fv["x"]:
            a = fv["x"][sid]["mcp"]
            if a["standing"] == "not exercised":
                classes[a["class_word"]] = classes.get(a["class_word"], 0) + 1
        text = "%s %d" % (fv["title"], n)
        if len(classes) > 1:
            text += " (%s)" % ", ".join("%d %s" % (c, w) for w, c in sorted(classes.items(), key=lambda kv: -kv[1]))
        parts.append(text)
    # the counts of the tests not exercised via MCP, red wherever printed (the owner's decision of 2026-10-16)
    v["mcp_not_exercised_by_family"] = cells.attn_marked(", ".join(parts))
    # the not-exercised tests split by why (the delta verifier's M2): an answer-side scenario-level mark (the request
    # expressible, the MCP answer lacking the field the property checks read: blocks.answer_mark_text) apart from
    # the requests the tools cannot express; `n_mcp_answer_lacks` is 0 for a manifest without such a mark
    lacks, tools_parts = [], []
    for fid in fams:
        fv = fam[fid]
        ids = [sid for sid in fv["x"] if fv["x"][sid]["mcp"]["standing"] == "not exercised"]
        mine = [sid for sid in ids if fv["x"][sid]["mcp"]["class"] == "answer_lacks_field"]
        lacks += mine
        if len(ids) > len(mine):
            tools_parts.append("%s %d" % (fv["title"], len(ids) - len(mine)))
    fields = []
    for sid in lacks:
        for _, check, _, _ in blocks.scenario_mark_parts(run.scenario(sid)):
            for f in blocks.mark_fields(run.scenario(sid), check) or [check]:
                if f not in fields:
                    fields.append(f)
    v["n_mcp_answer_lacks"] = len(lacks)
    v["mcp_answer_lacks_ids"] = ", ".join(lacks) or "none"
    v["mcp_answer_lacks_fields"] = ", ".join(fields) or "none"
    v["n_mcp_not_exercised_tools"] = sum(fam[fid]["n_mcp_not_exercised"] for fid in fams) - len(lacks)
    v["mcp_not_exercised_tools_by_family"] = cells.attn_marked(", ".join(tools_parts))
    n_exec, n_tests, text = _findings(ctx, sids)
    v["n_finding_executions"], v["n_finding_tests"], v["findings_text"] = n_exec, n_tests, text
    totals = run.summary().get("totals") or {}
    run_findings = sum((totals.get(a) or {}).get("outcomes", {}).get("finding", 0) for a in ARMS)
    v["run_findings_only_here"] = bool(n_exec) and n_exec == run_findings
    v["no_200_over_budget"] = not any(
        l.get("http") == 200 and l.get("outcome") == "finding" for sid in sids for arm in ARMS for l in run.arm_lines(sid, arm))
    v.update(_late503(ctx, sids))
    v.update(_stops(ctx, sids))
    v.update(_budget_stops(ctx, sids))
    # the timing, from the summary's own figures (the chapter's summary block reads the same)
    timing = (run.summary().get("timing") or {}).get(ch) or {}
    for arm in ARMS:
        k = ARM_KEY[arm]
        t = (timing.get(arm) or {}).get("by_host_cache") or {}
        for state in ("cold", "warm"):
            w = (t.get(state) or {}).get("wall") or {}
            s = (t.get(state) or {}).get("server") or {}
            v["wall_p50_%s_%s" % (state, k)] = cells.fmt_exec_text(w.get("p50")) if w.get("n") else "–"
            v["wall_p90_%s_%s" % (state, k)] = cells.fmt_exec_text(w.get("p90")) if w.get("n") else "–"
            v["n_%s_%s" % (state, k)] = w.get("n") or 0
            v["server_p50_%s_%s" % (state, k)] = cells.fmt_exec_text(s.get("p50")) if s.get("n") else "–"
            v["n_server_%s_%s" % (state, k)] = s.get("n") or 0
    backends = [l.get("mcp_phases", {}).get("backend_ms") for sid in sids for l in run.arm_lines(sid, "mcp")
                if isinstance(l.get("mcp_phases"), dict) and cells.host_cache_of(l) == "warm"
                and cells._num(l.get("mcp_phases", {}).get("backend_ms")) is not None]
    v["backend_p50_warm_mcp"] = cells.fmt_exec_text(cells.percentile(backends, 50)) if backends else "–"
    v["n_backend_warm_mcp"] = len(backends)
    # the chapter's executions in the host's memory episodes, for its Timing paragraph (RV1, V9; sections, read late:
    # it imports this module)
    import sections as _sections
    v["memory_note"] = _sections.memory_chapter_note(ctx, ch) if getattr(ctx, "run", None) is not None else ""
    # manifest v8: chapter G's five-way comparison and chapter S's stop records (doc/graphlets.py)
    if getattr(ctx, "run", None) is not None:
        import graphlets as _graphlets
        v.update(_graphlets.chapter_values(ctx, ch))
    cold_by_family = [(fam[fid]["n_cold_dl"], fid) for fid in fams if fam[fid]["n_cold_dl"]]
    v["cold_dl_families"] = ", ".join("%s %d" % (fid, n) for n, fid in sorted(cold_by_family, key=lambda t: (-t[0], t[1])))
    held, unsupported, items = _claims(ctx, ch, fam)
    v["claims_held"] = "; ".join(held)
    v["claims_unsupported"] = "; ".join(unsupported)
    v["n_claims_held"], v["n_claims_unsupported"] = len(held), len(unsupported)
    v["claims_items"] = items
    # the MCP agreements short of fully checked (READ-110), split by what they compared (the readability review's
    # verification of 2026-10-16, chapters issue 2): on the generic checks or a refusal alone (no property check
    # compared: the agreement verifies no property), or with some property check compared and agreeing on both
    # access paths and another not compared on some repeat (the detail block names which); the first such pair in
    # run order named with what it left uncompared, for the chapter summary's reading line
    v["n_mcp_agree_partial"] = max(0, v["n_mcp_agree"] - v["n_mcp_full"])
    v.update(_mcp_agree_short(run, sids, footnoted=lambda sid: bool(x[sid]["footnoted"])))
    v["mcp_agree_partial_example"] = ""
    for sid in sids:
        e = run.comp_entry(sid)
        d = (run.depth or {}).get(sid)
        if (e.get("mcp_agreement") or {}).get("class") == "agree" and d and d.get("class") != "full":
            nct = blocks.not_compared_text(d, list(run.check_marks(sid)) if hasattr(run, "check_marks") else [])
            v["mcp_agree_partial_example"] = "%s (%s%s)" % (
                sid, blocks.DEPTH_WORD.get(str(d.get("class")), str(d.get("class"))), ("; " + nct) if nct else "")
            break
    footnoted = [sid for sid in sids if x[sid]["footnoted"]]
    v["n_footnoted"] = len(footnoted)
    v["footnoted_ids"] = ", ".join(footnoted)
    v["footnoted_text"] = "; ".join(_members_phrase(x[sid]) for sid in footnoted)
    gap_fams = set((ctx.prose.get(ch) or {}).get("coverage_gap_families") or [])
    elsewhere = [sid for sid in footnoted if x[sid]["family"] not in gap_fams]
    v["n_footnoted_elsewhere"] = len(elsewhere)
    # the footnoted tests of the families the coverage gap discusses itself: when there are none, 'outside <the
    # family> every test passed' would imply a non-pass inside it (VERIFY-COLOURS-V5 item 7)
    v["n_footnoted_in_gap"] = len(footnoted) - len(elsewhere)
    v["footnoted_elsewhere_text"] = "; ".join(_members_phrase(x[sid]) for sid in elsewhere)
    # the same tests by id with their variants alone (READ-111: the chapter page points at the exception index
    # and the detail blocks instead of repeating each execution's chronology)
    v["footnoted_elsewhere_ids"] = ", ".join(_label(x[sid]) for sid in elsewhere)
    blocking = _gate_blocking(ctx, footnoted)
    v["n_footnoted_gate_blocking"] = len(blocking)
    v["footnoted_all_gate_blocking"] = bool(footnoted) and len(blocking) == len(footnoted)
    v["gate_judged"] = isinstance(run.gate, dict)
    v["logan_caveat"] = logan_caveat(run)
    # the chapter's attention counts (the contents' tag and I.1's attention table; the owner's decision of 2026-10-16)
    v.update(attention_facts(ctx, sids))
    return v


# -- needs attention (the owner's decision of 2026-10-16) ------------------------------------------------------------
# Red means "needs attention". A chapter's attention tag (the contents, under the chapter's line) and I.1's attention
# table count TESTS per kind, a test counted once per kind and possibly under several kinds, in this order: the
# defects first (a failure, an error, a finding, a disagreement between access paths), then the time-budget issues,
# then the tests not implemented via MCP. The kinds are defined once, in the reading guide's legend
# (sections.attention_legend_sentence, from ATTENTION_DEFINITION).
ATTENTION_KINDS = (("failure", "failure", "failures"), ("error", "error", "errors"), ("finding", "finding", "findings"),
                   ("disagreement", "disagreement", "disagreements"),
                   ("time_budget", "time-budget issue", "time-budget issues"),
                   ("not_implemented", "not implemented via MCP", "not implemented via MCP"))
ATTENTION_DEFINITION = {
    "failure": "a test with an execution that failed a check, on any access path",
    "error": "a test with an execution that got no HTTP answer or whose evaluation raised",
    "finding": "a test with an execution that broke a stated bound on the clock its check reads",
    "disagreement": "a test whose access paths disagree (the comparison's class disagree, direct or via MCP)",
    "time_budget": "a test with an execution that hit its time budget (a 503 deadline, an answer the deadline "
                   "limited or withheld, a walk or a resolve stopped by time, the attempt's duration bound, the "
                   "service's own deadline) or whose wall prints above %s",
    "not_implemented": "a test whose MCP access path stands not implemented via MCP (the tools cannot express the "
                       "request), or one of whose MCP executions had a property check set aside by the manifest's "
                       "property-level mark not implemented via MCP",
}
# the not-implemented definition where the run has an answer-side scenario-level mark (the delta verifier's M2)
ATTENTION_NOT_IMPLEMENTED_ANSWER_SIDE = (
    "a test whose MCP access path stands not implemented via MCP (the tools cannot express the request, or, by the "
    "manifest's scenario-level mark, the MCP answer lacks every field its property checks read), or one of whose MCP "
    "executions had a property check set aside by the manifest's property-level mark not implemented via MCP")
_SERVICE_DEADLINES = ("resolve_timeout", "traversal_resolve_deadline")


def time_budget_issue(line):
    """Whether an execution hit its time budget (ATTENTION_DEFINITION['time_budget']): a 503 deadline, a recorded
    time limitation (the runner's `time_limited`: an answer the deadline limited or withheld, a seed or a resolve
    stopped by time, the attempt's duration bound, the host's deadline relayed by the service), the service's own
    deadline read as a host state, or a wall printed above 55 s (red by the printed value, cells.time_tint)."""
    code, _ = _err(line)
    if line.get("http") == 503 and code == "deadline":
        return True
    if str(line.get("time_limited") or ""):
        return True
    if line.get("inconclusive_reason") in _SERVICE_DEADLINES:
        return True
    return cells.time_tint(line.get("wall_ms")) == "red"


def not_implemented_issue(run, sid, st):
    """Whether a test is not implemented via MCP (ATTENTION_DEFINITION['not_implemented']): its MCP arm standing ○,
    or an MCP execution whose checks the manifest's property-level mark not_implemented sets aside (read
    not_exposed; suite_common.checks_marked_aside), or that a scenario-level mark of not-implemented checks left
    inconclusive."""
    if st.present and st.glyph_key == "not_implemented_mcp":
        return True
    return cells.property_not_implemented(run, sid, st.lines)


def attention_of(ctx, sids):
    """{kind: [sids]} of the ATTENTION_KINDS over the scenarios (each test once per kind, in the given order)."""
    run = ctx.run
    out = dict((k, []) for k, _, _ in ATTENTION_KINDS)
    for sid in sids:
        states = ctx.state(sid)
        kinds = set()
        for arm in ARMS:
            st = states[arm]
            for l in st.lines:
                o = l.get("outcome")
                if o in ("fail", "other"):
                    kinds.add("failure")
                elif o == "error":
                    kinds.add("error")
                elif o == "finding":
                    kinds.add("finding")
                if time_budget_issue(l):
                    kinds.add("time_budget")
            if st.ne:
                kinds.add("disagreement")
        if not_implemented_issue(run, sid, states["mcp"]):
            kinds.add("not_implemented")
        for k in kinds:
            out[k].append(sid)
    return out


def attention_items(counts):
    """[(kind, n, words)] of the nonzero kinds in ATTENTION_KINDS' order: ('failure', 2, '2 failures')."""
    out = []
    for k, one, many in ATTENTION_KINDS:
        n = counts.get(k, 0)
        if n:
            out.append((k, n, "%d %s" % (n, one if n == 1 else many)))
    return out


def attention_text(counts, property_ids=None):
    """'2 time-budget issues, 9 not implemented via MCP': the nonzero kinds in order; '' when none. `property_ids`:
    the tests counted under not implemented via MCP by a property-level mark alone (their MCP access path does not
    stand ○), named after the count so that the reader finds them (VERIFY-MERGE L1): '18 not implemented via MCP
    (3 by a property-level mark: T-labels-0003, -0005 and -0006)'."""
    out = []
    for k, n, w in attention_items(counts):
        if k == "not_implemented" and property_ids:
            w += " (%d by a property-level mark: %s)" % (len(property_ids), _ids_short(property_ids))
        out.append(w)
    return ", ".join(out)


def attention_facts(ctx, sids):
    """The attention counts of a chapter (or of any set of tests) as facts: n_attention_<kind> (tests),
    attention_<kind>_ids, attention_text (the nonzero kinds in words, '' when none), n_attention_tests (the tests
    under any kind), attention_property_ids (the tests not implemented via MCP by a property-level mark alone, their
    MCP access path not standing ○; ArmState.ni_marked) and n_attention_property."""
    by = attention_of(ctx, sids)
    v = {}
    for k, _, _ in ATTENTION_KINDS:
        v["n_attention_%s" % k] = len(by[k])
        v["attention_%s_ids" % k] = ", ".join(by[k])
    prop = [sid for sid in by["not_implemented"] if ctx.state(sid)["mcp"].glyph_key != "not_implemented_mcp"]
    v["attention_property_ids"] = ", ".join(prop)
    v["n_attention_property"] = len(prop)
    v["attention_text"] = attention_text(dict((k, len(by[k])) for k in by), prop)
    v["n_attention_tests"] = len(set(sid for ids in by.values() for sid in ids))
    return v


AGREE_SHORT_GENERIC = ("generic_only", "refusal")
AGREE_SHORT_PARTIAL = ("partial", "partial_non_pass")
AGREE_SHORT_ON_FAILURE = ("on_failure",)


def _mcp_agree_short(run, sids, footnoted=None):
    """The agreeing MCP pairs short of fully checked, by depth class (report.DEPTH_CLASSES): {n_mcp_agree_generic
    (generic checks only or a refusal: no property check compared), n_mcp_agree_partial_checked (partial or
    partial_non_pass: some property check compared and agreeing on both, another not compared on some repeat),
    n_mcp_agree_on_failure (every property check compared, on a failing answer), mcp_agree_short_text (the three in
    words, separated by semicolons since the first holds a relative clause (the re-verification L-j), the empty ones
    left out; '' when none)}."""
    counts = {"generic": 0, "partial": 0, "on_failure": 0}
    partial_ids = []
    for sid in sids:
        e = run.comp_entry(sid)
        d = (run.depth or {}).get(sid)
        if (e.get("mcp_agreement") or {}).get("class") != "agree" or not d or d.get("class") == "full":
            continue
        cls = str(d.get("class"))
        if cls in AGREE_SHORT_GENERIC:
            counts["generic"] += 1
        elif cls in AGREE_SHORT_PARTIAL:
            counts["partial"] += 1
            partial_ids.append(sid)
        elif cls in AGREE_SHORT_ON_FAILURE:
            counts["on_failure"] += 1
    parts = []
    if counts["generic"]:
        parts.append("%d on the generic checks or a refusal alone, which verifies no property" % counts["generic"])
    if counts["partial"]:
        # the detail block names the uncompared check only where the test has one (`footnoted`: the chapter's
        # footnoted tests); else the tests and their checks are named here (the verification of round 2, R5)
        if footnoted is None or all(footnoted(sid) for sid in partial_ids):
            which = "the detail block names which"
        else:
            by_reason = []
            for sid in partial_ids:
                why = blocks.not_compared_text((run.depth or {}).get(sid) or {}, list(run.check_marks(sid))
                                               if hasattr(run, "check_marks") else []) or "a check not compared"
                for r in by_reason:
                    if r[0] == why:
                        r[1].append(sid)
                        break
                else:
                    by_reason.append([why, [sid]])
            which = "; ".join("%s: %s" % (ids_ranges(ids), why) for why, ids in by_reason)
        parts.append("%d with some property check compared and agreeing on both access paths and another not "
                     "compared on some repeat (%s)" % (counts["partial"], which))
    if counts["on_failure"]:
        parts.append("%d with every property check compared on a failing answer" % counts["on_failure"])
    return {"n_mcp_agree_generic": counts["generic"], "n_mcp_agree_partial_checked": counts["partial"],
            "n_mcp_agree_on_failure": counts["on_failure"], "mcp_agree_short_text": "; ".join(parts)}


def logan_caveat(run):
    import report
    lf = run.load or {}
    return bool((lf.get("sra_mean") or 0) > report.LOAD_CAVEAT_SRA_CPU_PCT
                or (lf.get("load_mean") or 0) > report.LOAD_CAVEAT_LOAD1)
