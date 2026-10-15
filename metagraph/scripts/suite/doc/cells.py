"""The cell grammar (DESIGN-DOC §3): how a scenario's standing on one arm becomes a glyph, its marks and its times;
what "clean" means; the aggregates of replicates cells and summary lines; the footnote numbering.

An `ArmState` is computed once per (scenario, arm) from the execution lines and compare's per-scenario entry; the
layouts in sections.py only format it.
"""

import os
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


def fmt_delta(ms):
    """A signed difference of two durations, rounded as fmt_t rounds: `−7 ms`, `+0.2 s`; below half a millisecond
    the sign and one decimal (`+0.2 ms`), so that a difference never loses its sign (I.3's rule); `0 ms` for an
    exact zero; `–` for none. Used wherever a difference is printed (a delta, an overhead, a net cost)."""
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


def fmt_times(values, marks=None):
    """Several durations in order, joined with ` / `, the unit once when they share it (`3.1 / 0.9 s`,
    `26 ms / 1.2 s`); `marks[i]` (`c`, `a`) follows the i-th number as a superscript."""
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
        if footnote:
            out += "\\fn{%d}" % footnote
        return out


def arm_state(run, sid, arm):
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
            st.reasons.append("the direct arms disagree")
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
    # the members the arm never exercised (the service refused or could not express the request, no execution
    # recorded) are said so, never counted as ok 0/0 (the owner's review, point 4)
    extras = []
    if ni:
        extras.append("%d\\,\\gni{}" % ni)
    if bd:
        extras.append("%d\\,\\gbd{}" % bd)
    if missing:
        extras.append("%d\\,\\gnone{}" % missing)
    not_exercised = bd + ni + missing
    denominator = n - not_exercised
    parts = ["ok %d/%d" % (ok, denominator)] if denominator else []
    if not_exercised:
        parts.append("not exercised: %s" % (extras[0] if len(extras) == 1 else "%d (%s)" % (
            not_exercised, " ".join(extras))))
    if walls:
        wall_text = fmt_pair(percentile(walls, 50), max(v for v in walls if _num(v) is not None))
        mark = host_marks_text(states)
        if mark:
            wall_text += "\\tmark{%s}" % mark
        parts.append(wall_text)
    text = " · ".join(parts) if parts else "–"
    srv = [v for v in servers if _num(v) is not None]
    line2 = "[%s]" % fmt_t(percentile(srv, 50)) if srv else ("[–]" if walls else None)
    return text, line2


def coverage_text(states, arm):
    """The coverage of a replicates unit on an arm, in words: `exercised 6/6, passed 6`, or `not exercised (3 ○)`
    when no member reached a tested execution, or `exercised 4/6, passed 3; not exercised 2 ○`."""
    members = [s[arm] for s in states.values()]
    n = len(members)
    exercised = [m for m in members if m.present and m.glyph_key not in ("service_by_design",
                                                                            "not_implemented_mcp")]
    passed = sum(1 for m in exercised if m.worst == "pass")
    rest = {}
    for m in members:
        if m in exercised:
            continue
        rest[m.glyph_key] = rest.get(m.glyph_key, 0) + 1
    rest_text = " ".join("%d\\,%s{}" % (rest[k], GLYPH.get(k, "\\gnone")) for k in GLYPH_ORDER if rest.get(k))
    if not exercised:
        return "not exercised (%s)" % rest_text if rest_text else "not exercised"
    text = "exercised %d/%d, passed %d" % (len(exercised), n, passed)
    if rest_text:
        text += "; not exercised %s" % rest_text
    return text


def property_tested(line, is_property=None):
    """What an execution tested of its scenario's property checks: ('yes' | 'partly' | 'no' | 'n/a', passed,
    total). A check is a property check by its recorded `property` flag, else by `is_property(name)`; a scenario
    without one (a refusal, whose status is the claim) reads n/a with its generic checks."""
    checks = [c for c in (line.get("checks") or []) if isinstance(c, dict)]
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
    "series": "times: wall per execution in repeat order, a / b (mcp: the cold call / the warm call); [ ] the "
              "host's server time; be the host call as the service timed it; c a i on a direct time: host-cold, "
              "aged, warm after an interrupted execution",
    "list": "times: wall per execution in repeat order, a / b (mcp: cold / warm); [ ] the host's server time; be "
            "the host call as the service timed it; c a i: host-cold, aged, warm after an interrupted execution",
    "pairs": "times: wall per execution in repeat order, a / b (mcp: cold / warm); [ ] the host's server time; "
             "c a i: host-cold, aged, warm after an interrupted execution; = and ≠ on a repeated test: equal / "
             "different answer digests",
    "replicates": "times: a / b is the median / the maximum of the members' walls (each member's first "
                  "host-warm execution, else its only one); [ ] the median server time",
    "refusals": "no times in the rows: the maximum wall per arm is in the summary line",
    "summary": "a / b in a summary line is p50 / p90",
}


def glyph_counts(states_list, arm):
    """`27 ● 1 ◐` for the states of several scenarios on an arm (LaTeX): each number bound to its glyph with a
    thin non-breaking space, a wider gap between the pairs; the glyph macro closed with {} so that the space
    after the last pair survives."""
    counts = {}
    for s in states_list:
        k = s[arm].glyph_key
        counts[k] = counts.get(k, 0) + 1
    parts = []
    for k in GLYPH_ORDER:
        if counts.get(k):
            parts.append("%d\\,%s{}" % (counts[k], GLYPH.get(k, "\\gnone")))
    return "\\enspace ".join(parts) if parts else "–"


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
        parts.append(fmt_t(max(vals)) if vals else "–")
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
