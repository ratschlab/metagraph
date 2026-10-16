"""The cell grammar (DESIGN-DOC §3): how a scenario's standing on one arm becomes a glyph, its marks and its times
(a time above a threshold on its tint: TIME_TINT_YELLOW_MS, TIME_TINT_RED_MS); what "clean" means; the aggregates
of replicates cells and summary lines; the footnote numbering.

An `ArmState` is computed once per (scenario, arm) from the execution lines and compare's per-scenario entry; the
layouts in sections.py only format it.
"""

import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
if os.path.dirname(HERE) not in sys.path:
    sys.path.insert(0, os.path.dirname(HERE))

import suite_common as sc  # noqa: E402
from latex import esc  # noqa: E402

ARMS = ("direct-local", "direct-remote", "mcp")
ARM_SHORT = {"direct-local": "dl", "direct-remote": "dr", "mcp": "mcp"}

# worst first
SEVERITY = {"error": 0, "fail": 1, "other": 2, "finding": 3, "inconclusive": 4, "mixed": 5, "service_by_design": 6,
            "not_implemented_mcp": 7, "pass": 8}
GLYPH = {"pass": "\\gpass", "fail": "\\gfail", "error": "\\gerr", "finding": "\\gfind", "inconclusive": "\\ginc",
         "not_implemented_mcp": "\\gni", "service_by_design": "\\gbd", "other": "\\gfail", None: "\\gnone"}
GLYPH_CHAR = {"pass": "●", "fail": "✖", "error": "⊘", "finding": "▲", "inconclusive": "◐", "not_implemented_mcp": "○",
              "service_by_design": "◆", "other": "✖", None: "–"}
WORD = {"pass": "pass", "fail": "fail", "error": "error (no HTTP answer)", "finding": "finding",
        "inconclusive": "inconclusive", "not_implemented_mcp": "not implemented via MCP",
        "service_by_design": "service-side by design", "other": "outcome outside the vocabulary", None: "no execution"}
GLYPH_ORDER = ("pass", "fail", "error", "finding", "inconclusive", "not_implemented_mcp", "service_by_design",
               "other", None)
WARM_AFTER_INTERRUPTED = "warm_after_interrupted"


# -- times -------------------------------------------------------------------------------------------------------

def _num(v):
    return v if isinstance(v, (int, float)) and not isinstance(v, bool) else None


def fmt_t(ms):
    """A duration (never negative): `26 ms` under 100 ms, `3.1 s` under 10 s, `18 s` above; `–` for none.
    The rounding, half up: whole milliseconds below 100 ms, one decimal in seconds below 10 s, whole seconds
    above (150 to 249 ms print as 0.2 s)."""
    v = _num(ms)
    if v is None:
        return "–"
    if v < 0.5:
        return "<1 ms"
    if v < 99.5:
        return "%d ms" % int(v + 0.5)
    if v < 9950:
        return "%.1f s" % (int(v / 100.0 + 0.5) / 10.0)
    return "%d s" % int(v / 1000.0 + 0.5)


# The tints of a printed execution time (the owner's decision of 2026-10-15): a wall, a server time, a backend time
# or an MCP phase, alone or as a median, a p90 or a maximum, strictly above TIME_TINT_YELLOW_MS is set on a yellow
# tint and strictly above TIME_TINT_RED_MS on a red one (doc.sty's timeyellow / timered), decided by the value as
# printed (the owner's decision of 2026-10-16: '10 s' is never tinted, '56 s' always red). The thresholds are this
# document's own: no budget, timing class or window of the suite is read for them, and nothing else reads them.
# A difference (fmt_delta: a delta, an overhead, how late a 503 was), a bound (a budget, a cap) and a distance in
# time (started after, ran before) carry none; a figure quoted in prose or inside a check's text carries its tint
# as a marker that latex.esc sets (fmt_exec_text, tint_quoted_times, esc_times).
TIME_TINT_YELLOW_MS = 10000.0
TIME_TINT_RED_MS = 55000.0
TINT_MACRO = {"yellow": "\\tyellow", "red": "\\tred"}


def printed_ms(text):
    """The duration a printed time shows, in milliseconds: '26 ms' 26, '3.1 s' 3100, '18 s' 18000, '15,593 ms'
    15593, '<1 ms' 0; None for '–' or a text that is no duration."""
    m = re.match(r"^\s*(<)?(\d[\d,]*(?:\.\d+)?)\s*(ms|s)\s*$", str(text or ""))
    if not m:
        return None
    if m.group(1):
        return 0.0
    v = float(m.group(2).replace(",", ""))
    return v * 1000.0 if m.group(3) == "s" else v


def tint_of(printed):
    """The tint of a printed duration's value in milliseconds (printed_ms): 'red' above TIME_TINT_RED_MS, 'yellow'
    above TIME_TINT_YELLOW_MS, else None."""
    v = _num(printed)
    if v is None:
        return None
    if v > TIME_TINT_RED_MS:
        return "red"
    if v > TIME_TINT_YELLOW_MS:
        return "yellow"
    return None


def time_tint(ms):
    """The tint of the duration `ms` as fmt_t prints it, decided by the printed value (the owner's decision of
    2026-10-16, what the reader sees; report.py marks its own printed times the same way, report.fmt_time, since
    the merge with gr/suite-v6): 10,400 ms prints '10 s' and has none, 10,600 ms
    prints '11 s' and is yellow, 55,400 ms prints '55 s' and is yellow; None for no value."""
    if _num(ms) is None:
        return None
    return tint_of(printed_ms(fmt_t(ms)))


def tint_marked(text, printed=None):
    """`text` (plain, the printed form of an execution time) carrying its tint as a marker that latex.esc sets,
    the tint decided by the printed value (`printed`, in ms, else read from `text`); `text` itself when none."""
    import latex
    tint = tint_of(printed if printed is not None else printed_ms(text))
    return "%s%s%s" % (latex.TINT_OPEN[tint], text, latex.TINT_CLOSE) if tint else text


def tinted(text, ms):
    """`text` (LaTeX, the printed form of the duration `ms`) on its tint when it has one: `\\tred{61 s}`."""
    tint = time_tint(ms)
    return "%s{%s}" % (TINT_MACRO[tint], text) if tint else text


def fmt_exec(ms):
    """An execution time as fmt_t prints it, on its tint when it has one (the cells, the detail blocks' rows, the
    aggregates and the timing tables print their times through this, never through fmt_t directly)."""
    return tinted(fmt_t(ms), ms)


# A figure quoted in prose or inside a check's text is plain text that latex.esc escapes afterwards, so its tint
# travels as a marker (latex.TINT_OPEN / TINT_CLOSE) that esc turns into the macro last. The execution times a
# check's or the runner's text quotes: budget_respected's 503 evidence (`client wall 10452 ms`, `a client wall time of
# 24281 ms`) and the host's own elapsed time it reads (`elapsed 388 ms`), and the MCP runner's host-state rule
# (`resolve_timeout after 30395 ms`, `traversal_resolve_deadline after 29201 ms`); a lateness
# (`delivered 7452 ms after`), a budget (`the 3000 ms budget`) and every other figure in such a text stay plain.
_QUOTED_TIME_RE = re.compile(r"(client wall (?:time of )?|(?:deadline|timeout) after |\belapsed )(\d[\d,]*(?:\.\d+)?) ms")


def fmt_exec_text(ms):
    """An execution time as fmt_t prints it, carrying its tint as a marker for a figure quoted in prose (latex.esc
    sets the macro): the facts the prose files fill in print their times through this."""
    import latex
    tint = time_tint(ms)
    text = fmt_t(ms)
    return "%s%s%s" % (latex.TINT_OPEN[tint], text, latex.TINT_CLOSE) if tint else text


def tint_quoted_times(text):
    """`text` (plain) with every execution time it quotes as a check's or the runner's text does (_QUOTED_TIME_RE)
    carrying its tint as a marker; the other figures of the text stay."""
    import latex

    def mark(m):
        try:
            ms = float(m.group(2).replace(",", ""))
        except ValueError:
            return m.group(0)
        tint = tint_of(ms)  # the figure as the text prints it is the printed value
        if not tint:
            return m.group(0)
        return "%s%s%s ms%s" % (m.group(1), latex.TINT_OPEN[tint], m.group(2), latex.TINT_CLOSE)
    return _QUOTED_TIME_RE.sub(mark, str(text or ""))


def esc_times(text, limit=None):
    """latex.esc of a check's or the runner's text with the execution times it quotes on their tints
    (tint_quoted_times): the detail blocks' finding, checks and host-state lines print such texts through this."""
    import latex
    return latex.esc(tint_quoted_times(text), limit)


# -- needs attention -----------------------------------------------------------------------------------------------
# Red means "needs attention" (the owner's decision of 2026-10-16). A test not implemented via MCP (the MCP-only
# class: the tools cannot express the request) is marked red wherever it is printed: the ○ glyph (doc.sty's \gni),
# its counts in the standing, outcomes, agreement and coverage lines, the property-level marks and the class words
# of the detail blocks; the attention tags of the contents and I.1's attention table are red too.

def attn(tex):
    """`tex` (LaTeX) set in the attention red (doc.sty's \\attn)."""
    return "\\attn{%s}" % tex if tex else tex


def attn_marked(text):
    """`text` (plain) carrying the attention marker that latex.esc turns into \\attn{…} (a phrase quoted in prose)."""
    import latex
    return "%s%s%s" % (latex.ATTN_OPEN, text, latex.TINT_CLOSE) if text else text


def fmt_delta(ms):
    """A signed difference of two durations, rounded as fmt_t rounds: `−7 ms`, `+0.2 s`; below half a millisecond
    the sign and one decimal (`+0.2 ms`), so that a difference never loses its sign (I.3's rule); `0 ms` for an
    exact zero; `–` for none. Used wherever a difference is printed (a delta, an overhead, a net cost); a
    difference is never tinted."""
    v = _num(ms)
    if v is None:
        return "–"
    if v == 0:
        return "0 ms"
    sign = "−" if v < 0 else "+"
    if abs(v) < 0.5:
        return "%s%.1f ms" % (sign, abs(v))
    return sign + fmt_t(abs(v))


def _split_unit(text):
    if text == "–" or text == "<1 ms":
        return text, None
    num, unit = text.rsplit(" ", 1)
    return num, unit


def fmt_times(values, marks=None, tint=True):
    """Several durations in order, joined with ` / `, the unit once when they share it (`3.1 / 0.9 s`,
    `26 ms / 1.2 s`); `marks[i]` (`c`, `a`) follows the i-th number as a superscript. Each figure sits on its own
    tint when its value has one (`\\tred{61} / 0.9 s`; the mark outside the tint); `tint=False` for a caller that
    prints something other than execution times."""
    values = list(values or [])
    if not values:
        return "–"
    texts = [fmt_t(v) for v in values]
    marks = list(marks or [None] * len(values))
    units = set(_split_unit(t)[1] for t in texts)
    parts = []
    for i, t in enumerate(texts):
        num, unit = _split_unit(t)
        shown = num if (len(units) == 1 and unit and i < len(texts) - 1) else t
        if tint:
            shown = tinted(shown, values[i])
        if marks[i]:
            shown += "\\tmark{%s}" % marks[i]
        parts.append(shown)
    return " / ".join(parts)


def fmt_pair(a, b):
    """`p50 / p90` with the unit once."""
    return fmt_times([a, b])


def fmt_n(v, digits=0):
    x = _num(v)
    if x is None:
        return "–"
    return ("%%.%df" % digits) % x


def fmt_gib(v):
    x = _num(v)
    return "–" if x is None else "%.0f GiB" % (x / float(1 << 30))


def percentile(values, q):
    return sc.percentile([v for v in values if _num(v) is not None], q)


# -- the state of a scenario on an arm --------------------------------------------------------------------------

def worst_of(outcomes):
    items = [o for o in outcomes if o is not None]
    if not items:
        return None
    return min(items, key=lambda o: SEVERITY.get(o, 2))


def host_cache_of(line):
    """cold / warm / warm_after_interrupted / aged / not_reached / unknown."""
    hc = line.get("host_cache") or line.get("cache") or "unknown"
    if hc == "warm" and line.get("host_cache_after_interrupted"):
        return WARM_AFTER_INTERRUPTED
    return hc


HOST_MARK = {"cold": "c", "aged": "a", WARM_AFTER_INTERRUPTED: "i"}


def host_mark(line):
    """The mark on a direct time: c host-cold, a aged, i host-warm after an interrupted execution; none when
    host-warm (or the host not reached, whose cell carries no time)."""
    return HOST_MARK.get(host_cache_of(line))


def host_marks_text(states):
    """The mark on an aggregate over executions in the host-cache states `states`: none when every one is
    host-warm, the letter alone when every one shares the state (`c`), else the letter with its share
    (`c 2/6`)."""
    n = len(states)
    if not n:
        return None
    parts = []
    for state, letter in (("cold", "c"), ("aged", "a"), (WARM_AFTER_INTERRUPTED, "i")):
        k = sum(1 for st in states if st == state)
        if k == n:
            parts.append(letter)
        elif k:
            parts.append("%s %d/%d" % (letter, k, n))
    return " ".join(parts) or None


class ArmState(object):
    def __init__(self, arm):
        self.arm = arm
        self.lines = []
        self.present = False
        self.missing_planned = False
        self.outcomes = []
        self.worst = None
        self.glyph_key = None
        self.glyph_marks = []  # m, t, s, n, u
        self.ne = False  # the disagreement mark
        self.clean = True
        self.reasons = []
        self.agreement = {}
        self.depth = None
        # an MCP arm not standing ○ one of whose executions had a property check set aside by the manifest's
        # property-level mark not implemented via MCP (property_not_implemented): the cell carries the small red ○
        # (doc.sty's \gnip) and the test counts under not implemented via MCP (facts.not_implemented_issue)
        self.ni_marked = False

    @property
    def glyph(self):
        return GLYPH.get(self.glyph_key, "\\gnone")

    @property
    def walls(self):
        return [l.get("wall_ms") for l in self.lines]

    @property
    def servers(self):
        return [((l.get("server") or {}).get("elapsed_ms") if isinstance(l.get("server"), dict) else None)
                for l in self.lines]

    @property
    def backends(self):
        return [((l.get("mcp_phases") or {}).get("backend_ms") if isinstance(l.get("mcp_phases"), dict) else None)
                for l in self.lines]

    @property
    def time_marks(self):
        if self.arm == "mcp":
            return [None] * len(self.lines)
        return [host_mark(l) for l in self.lines]

    def first_warm(self):
        """The first host-warm execution (else the first), for the aggregates."""
        for l in self.lines:
            if host_cache_of(l) == "warm":
                return l
        return self.lines[0] if self.lines else None

    def glyph_tex(self, footnote=None):
        out = self.glyph
        for m in self.glyph_marks:
            out += "\\tmark{%s}" % m
        if self.ne:
            out += "\\gne"
        if self.ni_marked:
            out += "\\gnip"
        if footnote:
            out += "\\fn{%d}" % footnote
        return out


def property_not_implemented(run, sid, lines):
    """Whether an MCP execution of the scenario among `lines` had a property check set aside by the manifest's
    property-level mark not implemented via MCP (the check read not_exposed: suite_common.checks_marked_aside, the
    reader's marked_aside), or was left inconclusive by a scenario-level mark over marks that include a
    not-implemented one (suite_common.SCENARIO_LEVEL_MARK). False for a manifest without property-level marks (v1 to
    v4). The one rule the cells' mark (ArmState.ni_marked) and the attention count (facts.not_implemented_issue)
    apply."""
    marks = run.check_marks(sid) if hasattr(run, "check_marks") else {}
    if not marks or not any(m.get("mark") == "not_implemented" for m in marks.values()):
        return False
    for l in lines:
        aside = run.marked_aside(sid, l) if hasattr(run, "marked_aside") else {}
        if any(m.get("mark") == "not_implemented" for m in aside.values()):
            return True
        if l.get("outcome") == "inconclusive" and l.get("inconclusive_reason") == sc.SCENARIO_LEVEL_MARK:
            return True
    return False


def arm_state(run, sid, arm):
    """The ArmState of a scenario on an arm (_arm_state), with the property-level mark of an MCP arm that does not
    stand ○ (ArmState.ni_marked)."""
    st = _arm_state(run, sid, arm)
    if arm == "mcp" and st.present and st.glyph_key != "not_implemented_mcp":
        st.ni_marked = property_not_implemented(run, sid, st.lines)
    return st


def _arm_state(run, sid, arm):
    st = ArmState(arm)
    lines = run.arm_lines(sid, arm)
    entry = run.comp_entry(sid)
    carm = run.comp_arm(sid, arm)
    da = entry.get("direct_agreement") or {}
    ma = entry.get("mcp_agreement") or {}
    st.agreement = ma if arm == "mcp" else da
    st.lines = lines
    if not lines:
        st.present = False
        st.missing_planned = run.planned_missing(sid, arm) or (arm == "mcp" and ma.get("class") == "missing")
        if st.missing_planned:
            st.clean = False
            st.reasons.append("a planned execution was not recorded")
        return st
    st.present = True
    st.outcomes = list(carm.get("outcomes") or [l.get("outcome") for l in lines])
    st.worst = worst_of(st.outcomes)
    if len(set(st.outcomes)) > 1:
        st.glyph_marks.append("m")
        st.reasons.append("the repeats' outcomes differ: %s" % " / ".join(str(o) for o in st.outcomes))
    st.glyph_key = st.worst if st.worst in GLYPH else "other"
    if arm == "mcp":
        cls = ma.get("class")
        if cls == "not_implemented_mcp" or st.worst == "not_implemented_mcp":
            st.glyph_key = "not_implemented_mcp"
            st.glyph_marks = [m for m in st.glyph_marks if m != "m"]
            return st
        if cls == "service_by_design" or (cls is None and st.worst == "service_by_design"):
            st.glyph_key = "service_by_design"
            st.glyph_marks = [m for m in st.glyph_marks if m != "m"]
            if ma.get("verified") is False:
                st.glyph_marks.append("u")
            return st
        if cls == "disagree":
            st.ne = True
            st.clean = False
            st.reasons.append("MCP disagrees with the reference arm")
        elif cls == "not_comparable":
            st.glyph_marks.append("n")
            st.clean = False
            st.reasons.append("not comparable: no check both arms evaluated")
        st.depth = (run.depth or {}).get(sid)
    elif arm == "direct-remote":
        cls = da.get("class")
        if cls == "disagree":
            st.ne = True
            st.clean = False
            st.reasons.append("the direct access paths disagree")
        elif cls in ("timing", "same_checks", "not_comparable"):
            st.glyph_marks.append({"timing": "t", "same_checks": "s", "not_comparable": "n"}[cls])
            st.clean = False
            st.reasons.append({"timing": "the arms differ by timing alone",
                               "same_checks": "different answers with the same check outcomes",
                               "not_comparable": "not comparable: no check both arms evaluated"}[cls])
    if st.worst != "pass":
        st.clean = False
        st.reasons.insert(0, WORD.get(st.worst, str(st.worst)))
    return st


def scenario_states(run, sid):
    return dict((arm, arm_state(run, sid, arm)) for arm in ARMS)


def is_clean(states):
    return all(s.clean for s in states.values())


SERVICE_DECIDED = ("service_by_design", "not_implemented_mcp")


def by_design_line(st, line):
    """Whether an execution of an MCP arm standing not exercised or by design is one the service decided in the
    host's place: every line of a not-exercised arm (a refusal the table files under not implemented reads not
    exercised whatever the runner wrote), and of a by-design arm the lines whose recorded outcome is the service's
    (service_by_design: the case's repeats, compare's `repeats`). The other repeat of a by-design arm (compare's
    `other_repeats`: a separate host answer the case did not match, a 200 that passed or was cut by time) keeps its
    recorded outcome (VERIFY-PLXM H2, VERIFY-QCND H2)."""
    if st.glyph_key == "not_implemented_mcp":
        return True
    if st.glyph_key == "service_by_design":
        return line.get("outcome") in SERVICE_DECIDED
    return False


def exercised_lines(st):
    """The executions of an arm the service did not decide in the host's place (by_design_line)."""
    return [l for l in st.lines if not by_design_line(st, l)]


# -- cells -------------------------------------------------------------------------------------------------------

def server_text(st):
    if not st.present:
        return None
    vals = st.servers
    if all(_num(v) is None for v in vals):
        return "[–]"
    return "[%s]" % fmt_times(vals)


def backend_text(st):
    if st.arm != "mcp" or not st.present:
        return None
    vals = st.backends
    if all(_num(v) is None for v in vals):
        return None
    return "be %s" % fmt_times(vals)


def cell_parts(st, footnote=None):
    """(line 1, line 2, line 3) of a cell with times; the by-design and not-implemented MCP cells carry the glyph
    alone (the time is the service's refusal, said in the appendix)."""
    if not st.present:
        return (st.glyph_tex(footnote), None, None)
    if st.glyph_key in ("not_implemented_mcp", "service_by_design"):
        return (st.glyph_tex(footnote), None, None)
    line1 = "%s\\ %s" % (st.glyph_tex(footnote), fmt_times(st.walls, st.time_marks))
    return (line1, server_text(st), backend_text(st))


def aggregate(states, arm):
    """The replicates cell of `states` ({sid: {arm: ArmState}}) on `arm`: ('ok n/n · p50 / max', '[server p50]')."""
    members = [s[arm] for s in states.values()]
    n = len(members)
    ok = sum(1 for m in members if m.present and m.worst == "pass" and m.clean)
    bd = sum(1 for m in members if m.glyph_key == "service_by_design")
    ni = sum(1 for m in members if m.glyph_key == "not_implemented_mcp")
    missing = sum(1 for m in members if not m.present)
    walls, servers, states = [], [], []
    for m in members:
        if not m.present or m.glyph_key in ("service_by_design", "not_implemented_mcp"):
            continue
        l = m.first_warm()
        if l is None:
            continue
        walls.append(l.get("wall_ms"))
        srv = (l.get("server") or {}).get("elapsed_ms") if isinstance(l.get("server"), dict) else None
        servers.append(srv)
        if arm != "mcp" and _num(l.get("wall_ms")) is not None:
            states.append(host_cache_of(l))
    # the members the arm never exercised (the service could not express the request, or no execution was
    # recorded) are said so, never counted as ok 0/0 (the owner's review, point 4); a member the service decided
    # by design is said apart (READ-208): the host may have been reached (a relayed refusal, a normalised
    # request) or not (a refusal at submit), which `by_design_reached_text` states
    extras = []
    if ni:
        extras.append(attn("%d\\,\\gni{}" % ni))
    if missing:
        extras.append("%d\\,\\gnone{}" % missing)
    not_exercised = bd + ni + missing
    denominator = n - not_exercised
    parts = ["ok %d/%d" % (ok, denominator)] if denominator else []
    if ni and not missing:
        parts.append(attn("not exercised: %d\\,\\gni{}" % ni))
    elif ni or missing:
        parts.append("not exercised: %s" % " ".join(extras))
    if bd:
        parts.append("by design: %d\\,\\gbd{}%s" % (bd, by_design_reached_text(
            [m for m in members if m.glyph_key == "service_by_design"])))
    marked = sum(1 for m in members if getattr(m, "ni_marked", False))
    if marked:
        parts.append(attn(PROPERTY_NI_COUNT % marked))
    if walls:
        wall_text = fmt_pair(percentile(walls, 50), max(v for v in walls if _num(v) is not None))
        mark = host_marks_text(states)
        if mark:
            wall_text += "\\tmark{%s}" % mark
        parts.append(wall_text)
    text = " · ".join(parts) if parts else "–"
    srv = [v for v in servers if _num(v) is not None]
    line2 = "[%s]" % fmt_exec(percentile(srv, 50)) if srv else ("[–]" if walls else None)
    return text, line2


def host_reached(line):
    """Whether an execution reached the host (its recorded host-cache state is one of the host's)."""
    return host_cache_of(line) in ("cold", "warm", WARM_AFTER_INTERRUPTED, "aged")


def by_design_reached_text(bd_members):
    """How many of the by-design members' service-decided executions reached the host, in a parenthesis for the
    grouped cell and the coverage line (READ-208): a relayed host refusal or a normalised request reached the
    host (the host executed the request; the service decided the pair), a refusal at submit did not. Empty when
    no member is by design."""
    if not bd_members:
        return ""
    reached = sum(1 for m in bd_members if any(host_reached(l) for l in m.lines if by_design_line(m, l)))
    if reached == 0:
        return " (decided before the host)"
    if reached == len(bd_members):
        return " (the host reached, the pair decided by the service)"
    return " (%d the host reached, the rest decided before it)" % reached


def coverage_text(states, arm):
    """The coverage of a replicates unit on an arm, in words: `exercised 6/6, passed 6`, or `exercised 0/3; not
    exercised 3 ○` when no member reached a tested execution, or `exercised 4/6, passed 3; not exercised 2 ○`; a
    member the service decided by design is said apart from one never exercised, with whether the host was
    reached (READ-208): `exercised 1/2, passed 1; by design 1 ◆ (the host reached, the pair decided by the
    service)`. A member standing by design whose other repeat the case did not match (by_design_line: a separate
    host answer) is exercised on that repeat: `exercised 14/22 (3 on one repeat only, the other by design),
    passed 9` (VERIFY-PLXM H2)."""
    members = [s[arm] for s in states.values()]
    n = len(members)
    exercised = [m for m in members if m.present and m.glyph_key not in SERVICE_DECIDED]
    partly = [m for m in members if m.present and m.glyph_key in SERVICE_DECIDED and exercised_lines(m)]
    passed = sum(1 for m in exercised if m.worst == "pass")
    passed += sum(1 for m in partly if all(l.get("outcome") == "pass" for l in exercised_lines(m)))
    rest = [m for m in members if m not in exercised and m not in partly]
    ni = [m for m in rest if m.present and m.glyph_key == "not_implemented_mcp"]
    bd = [m for m in rest if m.present and m.glyph_key == "service_by_design"]
    other = {}
    for m in rest:
        if m in ni or m in bd:
            continue
        other[m.glyph_key] = other.get(m.glyph_key, 0) + 1
    text = "exercised %d/%d" % (len(exercised) + len(partly), n)
    if partly:
        text += " (%d on one repeat only, the other by design)" % len(partly)
    if exercised or partly:
        text += ", passed %d" % passed
    if ni:
        text += "; " + attn("not exercised %d\\,\\gni{}" % len(ni))
    if bd:
        text += "; by design %d\\,\\gbd{}%s" % (len(bd), by_design_reached_text(bd))
    other_text = " ".join("%d\\,%s{}" % (other[k], GLYPH.get(k, "\\gnone")) for k in GLYPH_ORDER if other.get(k))
    if other_text:
        text += "; not recorded %s" % other_text
    marked = sum(1 for m in members if getattr(m, "ni_marked", False))
    if marked:
        text += "; " + attn(PROPERTY_NI_COUNT % marked)
    return text


_NOT_APPLICABLE_RE = re.compile(r"^http \S+: the fields are checked on an http \S+ only$")


def check_not_applicable(check):
    """Whether a recorded check item did not apply to the answer it was read on (checks.py's `answer_fields` with
    `when_http` on another status passes with the detail 'http 200: the fields are checked on an http 409 only'):
    nothing of it was tested, so it is neither counted as tested nor as untested (VERIFY-T M9)."""
    return bool(_NOT_APPLICABLE_RE.match(str(check.get("detail") or "")))


def property_tested(line, is_property=None, aside=None):
    """What an execution tested of its scenario's property checks: ('yes' | 'partly' | 'no' | 'n/a', passed,
    total). A check is a property check by its recorded `property` flag, else by `is_property(name)`; a scenario
    without one (a refusal, whose status is the claim) reads n/a with its generic checks. |aside|: the keys
    (check name, occurrence) of the items a property-level mark of the manifest sets aside on this MCP execution
    (suite_common.checks_marked_aside: marked not implemented via MCP or service-side by design, read not_exposed):
    they are neither tested nor untested here and are left out of the count (the caller names them)."""
    checks = [c for c in (line.get("checks") or []) if isinstance(c, dict) and not check_not_applicable(c)]
    if aside:
        seen = {}
        kept = []
        for c in checks:
            name = c.get("check")
            key = (name, seen.get(name, 0))
            seen[name] = seen.get(name, 0) + 1
            if key not in aside:
                kept.append(c)
        checks = kept
    props = [c for c in checks if (c.get("property") is True or (c.get("property") is None and is_property
                                                                   and is_property(c.get("check"))))]
    if not props:
        generic_passed = sum(1 for c in checks if c.get("outcome") == "pass")
        return ("n/a", generic_passed, len(checks))
    passed = sum(1 for c in props if c.get("outcome") == "pass")
    if passed == len(props):
        return ("yes", passed, len(props))
    if passed:
        return ("partly", passed, len(props))
    return ("no", 0, len(props))


# the local timing legends: what `/` and the brackets mean under each table type (the owner's review, point 6)
TIMING_LEGEND = {
    "series": "times: the wall time of each execution in repeat order, a / b = repeat 0 / repeat 1 (via MCP the cold "
              "call / the warm call); [ ] the host's own server time; be the host call as the service timed it; c a i "
              "on a direct time: host-cold, aged, warm after an interrupted execution",
    "list": "times: the wall time of each execution in repeat order, a / b = repeat 0 / repeat 1 (via MCP cold / "
            "warm); [ ] the host's own server time; be the host call as the service timed it; c a i: host-cold, "
            "aged, warm after an interrupted execution",
    "pairs": "times: the wall time of each execution in repeat order, a / b = repeat 0 / repeat 1 (via MCP cold / "
             "warm); [ ] the host's own server time; c a i: host-cold, aged, warm after an interrupted execution; = "
             "and ≠ on a repeated test: equal / different answer digests",
    "replicates": "times: a / b = the median / the maximum of the group's walls, one wall per test (each test's "
                  "first host-warm execution, else its only one), never cold / warm; n counts the group's tests; "
                  "[ ] the median server time",
    "refusals": "no times in the rows: the maximum wall per access path is in the summary line",
    "summary": "a / b in a summary line is p50 / p90: the median and the 90th percentile over the host-cold or "
               "host-warm executions with a wall time, whatever their outcome (an execution not reached, aged or warm "
               "after an interrupted one is left out)",
}
# the reader-facing key beside a unit's heading (READ-105): what the two numbers of a time cell mean in this layout
LAYOUT_KEY = {
    "series": "series · times: wall per repeat, r0 / r1 (mcp: cold / warm)",
    "list": "list · times: wall per repeat, r0 / r1 (mcp: cold / warm)",
    "pairs": "pairs · times: wall per repeat, r0 / r1 (mcp: cold / warm)",
    "replicates": "grouped · times: wall median / maximum over the group's tests; n = tests",
    "refusals": "refusals · no times in the rows",
}


def repeats_texts(n, mcp_cold_first=True):
    """(the heading key, the legend) of a pairs table whose tests are repeated `n` times (chapter D; READ-321):
    every repeat's wall is printed in execution order, so the key and the legend name the repeats r0 / r1 / ...,
    say which MCP call was host-cold as the record has it (`mcp_cold_first`: repeat 0 cold, the later ones warm),
    define the bracketed server time and the `be` line per execution, and tie = / ≠ to the normalised answer
    digests of the repeats on that access path."""
    reps = " / ".join("r%d" % i for i in range(n))
    later = "r1" if n == 2 else ("r1 and r2" if n == 3 else "r1 to r%d" % (n - 1))
    mcp = ("via MCP r0 is the host-cold call, %s the warm ones" % later) if mcp_cold_first else \
        "via MCP the detail blocks state each call's host state"
    key = "pairs · times: wall per repeat, %s (mcp: %s)" % (reps, ("r0 cold, %s warm" % later) if mcp_cold_first
                                                             else "see the detail blocks")
    legend = ("times: the wall time of each execution in execution order, %s (%s); [ ] the host's own server time "
              "per execution, in the same order; be the host call as the service timed it, per execution (MCP "
              "only); c a i on a direct time: host-cold, aged, warm after an interrupted execution; = and ≠: the "
              "normalised answer digests of the repeats on that access path are equal / differ (the first repeat "
              "is the baseline, a repeat a time budget limited is left out)" % (reps, mcp))
    return key, legend


def glyph_counts(states_list, arm):
    """`27 ● 1 ◐` for the states of several scenarios on an arm (LaTeX): each number bound to its glyph with a
    thin non-breaking space, a wider gap between the pairs; the glyph macro closed with {} so that the space
    after the last pair survives."""
    counts, marked = {}, {}
    for s in states_list:
        k = s[arm].glyph_key
        counts[k] = counts.get(k, 0) + 1
        if getattr(s[arm], "ni_marked", False):
            marked[k] = marked.get(k, 0) + 1
    parts = []
    for k in GLYPH_ORDER:
        if counts.get(k):
            pair = "%d\\,%s{}" % (counts[k], GLYPH.get(k, "\\gnone"))
            pair = attn(pair) if k == "not_implemented_mcp" else pair
            if marked.get(k):
                # the members of this glyph with a property check the manifest's property-level mark set aside as
                # not implemented via MCP (red, as every not-implemented count; the owner's decision of 2026-10-16)
                pair += "\\ " + attn(PROPERTY_NI_COUNT % marked[k])
            parts.append(pair)
    return "\\enspace ".join(parts) if parts else "–"


# the count of an arm's tests with a property check set aside by the manifest's property-level mark not implemented
# via MCP, after the glyph pair they stand under (glyph_counts, the aggregates, the coverage line)
PROPERTY_NI_COUNT = "(%d with a property check not implemented via MCP)"


def wall_buckets(run, sids, arm):
    """{cold: [wall], warm: [wall]} over the execution lines of `sids` on `arm` (warm after a cut and the host not
    reached excluded)."""
    out = {"cold": [], "warm": []}
    for sid in sids:
        for l in run.arm_lines(sid, arm):
            hc = host_cache_of(l)
            if hc in out and _num(l.get("wall_ms")) is not None:
                out[hc].append(l["wall_ms"])
    return out


def timing_text(run, sids):
    """`wall p50 / p90 host-cold: a · b · c; host-warm: a · b · c` over the arms."""
    cold, warm = [], []
    for arm in ARMS:
        b = wall_buckets(run, sids, arm)
        cold.append(fmt_pair(percentile(b["cold"], 50), percentile(b["cold"], 90)) if b["cold"] else "–")
        warm.append(fmt_pair(percentile(b["warm"], 50), percentile(b["warm"], 90)) if b["warm"] else "–")
    return "wall p50 / p90 host-cold: %s; host-warm: %s" % (" · ".join(cold), " · ".join(warm))


def max_wall_text(run, sids):
    parts = []
    for arm in ARMS:
        vals = [l.get("wall_ms") for sid in sids for l in run.arm_lines(sid, arm) if _num(l.get("wall_ms")) is not None]
        parts.append(fmt_exec(max(vals)) if vals else "–")
    return "wall max: " + " · ".join(parts)


class Footnotes(object):
    """Footnote marks of one chapter: superscript digits in reading order, one per test that is not clean."""

    def __init__(self):
        self.numbers = {}
        self.order = []

    def assign(self, sid):
        if sid not in self.numbers:
            self.numbers[sid] = len(self.numbers) + 1
            self.order.append(sid)
        return self.numbers[sid]

    def get(self, sid):
        return self.numbers.get(sid)


def esc_t(text):
    return esc(text)
