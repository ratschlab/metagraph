"""The suite's checks (DESIGN §8.2): pure functions of a scenario and a normalised answer, stdlib only, shared by
every arm so that a disagreement between arms is never a difference of checking code.

    evaluate(scenario, answer, arm, previous=None)
        -> [{"check": name, "outcome": "pass" | "fail" | "not_exposed" | "inconclusive" | "finding", "detail": str}]

`answer` is the NORMALISED answer: {"http": int, "body": obj | None, "error": {"code", "text"} | None} and, where
the arm has them, "headers" (lower-case names). `normalise_direct` builds it from a direct HTTP answer; the MCP
arm's normaliser maps its tool results into the same shape for the fields both expose. A check that needs a field
the MCP arm does not expose answers "not_exposed" on that arm (and "fail" on a direct arm, where the field is part
of the contract). `inconclusive`: the check could not test its property on this answer (an allowed non-200, a
withheld or incomplete list that hides the source, a stated resource stop before the property). A `finding` is a
stated bound or SPEC promise the host broke without failing the request (the time budget, budget_respected): the
report lists findings in their own section. scenario_outcome(results) gives the execution's outcome. `previous`
(optional) lists the `answer_digest`s of the earlier repeats of the same scenario and arm, for `deterministic`.

Every check named in a manifest's `expect` is a function registered in CHECKS; new checks are added here only.
Coordinates: truth["start"] is 0-based on the record's forward strand; host placements (`nt_coords`) are 1-based
inclusive; resolve runs and traverse intervals 0-based half-open (the SPECs' conventions).
"""

from __future__ import annotations

import gzip
import hashlib
import json
import math
import re
from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence, Tuple

K = 31
_COMP = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}
_NCBIEAA = {   # the tables the suite's peptides use (seqtools.NCBIEAA holds all of them)
    1: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    2: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
    3: "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    4: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    5: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
    6: "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    11: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    12: "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    27: "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    28: "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    31: "FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
}
_CODONS = [a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"]
_AMBIG = {"B": "DN", "Z": "EQ", "J": "IL"}
_BUDGET_STOP_REASONS = {"max_steps", "time"}
# the stop reasons that can leave a count inexact (SPEC-pattern-search §7.4: at_least "a stop interrupted
# discovery"; bounds after a stop in a deferred scan): the budgets of every phase (§7.6, §14.4, §19.8, §20.5) and the
# stop_at_threshold stops (§7.6, §12.1, §19.8, §20.9)
_INEXACT_STOP_REASONS = {"max_steps", "time", "max_annotation_work", "max_memory", "max_predicate_work",
                         "max_contexts", "max_anchors", "max_paths", "max_predicate_contexts"}
_ORIENT_RANK = {"+": 0, "forward": 0, "-": 1, "reverse": 1, "=": 2, "palindromic": 2}

PASS, FAIL, NOT_EXPOSED = "pass", "fail", "not_exposed"
# the host broke a stated bound or a SPEC promise without failing the request: reported in its own section of the
# report, never counted as a pass (budget_respected)
FINDING = "finding"
# the check could not test its property on this answer (an allowed refusal, a withheld or cut list that hides
# the source, an arm a stated resource limit stopped before the property): counted apart, never as a pass
INCONCLUSIVE = "inconclusive"
OUTCOMES = (PASS, FAIL, NOT_EXPOSED, INCONCLUSIVE, FINDING)


class Missing(Exception):
    """A field the check needs is absent from the answer."""


class Undecidable(Exception):
    """The answer states a reason the check relies on but not what would show it can have happened (a /pattern time
    stop on an answer that states no timing at all): the check reads inconclusive on every arm (evaluate), never a
    pass on an unjudged statement."""


# ---------------------------------------------------------------------------------------------- answers
def normalise_direct(http_status: int, body: Any, headers: Optional[Dict[str, str]] = None) -> dict:
    """The normalised answer of a direct HTTP call, as lane S2's runner makes it (suite_common.normalise_direct:
    {http, body, error, headers}). |body| is the body bytes (gzip bytes are recognised by their magic and
    undone), a str, or an already parsed object."""
    if isinstance(body, (bytes, bytearray)):
        raw = bytes(body)
        if raw[:2] == b"\x1f\x8b":
            raw = gzip.decompress(raw)
    elif isinstance(body, str):
        raw = body.encode("utf-8")
    else:
        raw = json.dumps(body).encode("utf-8") if body is not None else b""
    if _suite_common is not None:
        return _suite_common.normalise_direct(int(http_status), raw, headers)
    text = raw.decode("utf-8", "replace")
    parsed: Any = None
    if text.strip():
        try:
            parsed = json.loads(text)
        except ValueError:
            parsed = text
    error = None
    if not 200 <= int(http_status) < 300:
        if isinstance(parsed, dict):
            msg = parsed.get("error")
            if msg is not None and not isinstance(msg, str):
                msg = json.dumps(msg, sort_keys=True)
            error = {"code": parsed.get("code"), "text": msg}
        else:
            error = {"code": None, "text": text[:4000] if text else None}
    keep = {k.lower(): v for k, v in (headers or {}).items() if k.lower() in ("content-type", "content-encoding")}
    return {"http": int(http_status), "body": parsed, "error": error, "headers": keep}


try:   # lane S2's shared helpers: the runner digests with them, so the determinism check must too
    import suite_common as _suite_common
except ImportError:   # pragma: no cover - the vendored copy without it falls back to the same algorithm
    _suite_common = None

_CAPABILITY_VOLATILE = ("measured_", "observed_")


def _strip(o: Any, prefixes: tuple) -> Any:
    if isinstance(o, dict):
        return {k: _strip(v, prefixes) for k, v in o.items()
                if k != "timing" and not (prefixes and isinstance(k, str) and k.startswith(prefixes))}
    if isinstance(o, list):
        return [_strip(x, prefixes) for x in o]
    return o


def answer_digest(answer: dict, route: Optional[str] = None) -> str:
    """The answer digest of DESIGN §8.3, suite_common.answer_digest's: sha256 of the canonical JSON of the
    normalised answer without `headers` and every `timing` key (on a capabilities route also without the
    measured_* / observed_* keys)."""
    if _suite_common is not None:
        return _suite_common.answer_digest(answer, route)
    prefixes = ()
    if route and route.split(None, 1)[-1].split("?", 1)[0].rstrip("/").endswith("capabilities"):
        prefixes = _CAPABILITY_VOLATILE
    view = _strip({k: v for k, v in answer.items() if k != "headers"}, prefixes)
    text = json.dumps(view, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
    return hashlib.sha256(text).hexdigest()


# ---------------------------------------------------------------------------------------------- helpers
def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def _comp(b: str) -> str:
    return b.translate(_COMP)


def _get(obj: Any, path: str) -> Any:
    cur = obj
    for part in path.split("."):
        if isinstance(cur, dict):
            if part not in cur:
                raise Missing(path)
            cur = cur[part]
        elif isinstance(cur, list):
            try:
                cur = cur[int(part)]
            except (ValueError, IndexError):
                raise Missing(path)
        else:
            raise Missing(path)
    return cur


def _body(ans: dict) -> Any:
    b = ans.get("body")
    if b is None:
        raise Missing("body")
    return b


def _entry(ans: dict, i: int) -> dict:
    pats = _get(_body(ans), "patterns")
    if not isinstance(pats, list) or i >= len(pats):
        raise Missing("patterns[%d]" % i)
    return pats[i]


def _count(c: Any) -> Tuple[Optional[int], str]:
    if not isinstance(c, dict) or "relation" not in c:
        raise Missing("count")
    return c.get("value"), c["relation"]


def _request_pattern(sc: dict, i: int) -> Tuple[str, str]:
    p = sc["request"]["patterns"][i]
    for kind in ("dna", "iupac", "protein"):
        if kind in p:
            return kind, str(p[kind]).upper()
    raise Missing("request pattern")


def _pattern_length(kind: str, p: str) -> int:
    return 3 * len(p) if kind == "protein" else len(p)


def _primary_count_name(entry: dict, sc_req: dict) -> str:
    counts = entry.get("counts", {})
    if "selected" in counts and sc_req.get("predicate") is not None:
        return "selected"
    if "contexts" in counts:
        return "contexts"
    if sc_req.get("long_search") == "supported_paths" and "supported_paths" in counts:
        return "supported_paths"
    return "paths"


def source_contexts(t: dict, k: int = K) -> List[Tuple[str, int]]:
    """(kmer, offset) of every record k-mer that holds the truth window, within the truth's flank."""
    out = []
    L = t["length"]
    rel = t["start"] - t["flank_start"]
    for p in range(0, k - L + 1):
        a = rel - p
        if a < 0 or a + k > len(t["flank"]):
            continue
        out.append((t["flank"][a:a + k], p))
    return out


def _codon_sets(residue: str, table: int) -> set:
    aas = _NCBIEAA.get(table)
    if aas is None:
        raise Missing("genetic code %d not in checks" % table)
    m = dict(zip(_CODONS, aas))
    if residue == "X":
        return {c for c in _CODONS if m[c] != "*"}
    if residue in _AMBIG:
        return {c for c in _CODONS if m[c] in _AMBIG[residue]}
    return {c for c in _CODONS if m[c] == residue}


def _translate(dna: str, table: int) -> str:
    m = dict(zip(_CODONS, _NCBIEAA[table]))
    return "".join(m.get(dna[i:i + 3], "?") for i in range(0, len(dna) - 2, 3))


def expected_substitutions(kind: str, pattern: str, instance: str, strand: str, table: int = 1) -> Optional[List[dict]]:
    """The substitutions the host must state for |instance| against |pattern| (SPEC-pattern-search.md §21.8):
    positions (residues) in the coordinates of P; for '-' the instance is read in P's orientation. None when the
    instance is no instance of P at any level (§21.3, sense substitutes: a residue's substitute is a sense codon of the
    table coding another residue, so a stop codon never substitutes one, and X and * hold no substitution)."""
    if strand in ("-", "reverse"):
        read = revcomp(instance)
    else:
        read = instance
    out = []
    if kind == "protein":
        for r, aa in enumerate(pattern):
            codon = read[3 * r:3 * r + 3]
            if codon not in _codon_sets(aa, table):
                found = dict(zip(_CODONS, _NCBIEAA[table])).get(codon, "?")
                if aa in ("X", "*") or found == "*":
                    return None
                out.append({"residue": r, "expected": aa, "found": found, "codon": codon})
        return out
    for i, code in enumerate(pattern):
        if read[i] not in IUPAC[code]:
            out.append({"position": i, "expected": code, "found": read[i]})
    return out


def _instance_matches(kind: str, pattern: str, instance: str, strand: str, table: int) -> bool:
    read = revcomp(instance) if strand in ("-", "reverse") else instance
    if kind == "protein":
        return all(read[3 * r:3 * r + 3] in _codon_sets(aa, table) for r, aa in enumerate(pattern))
    return len(read) == len(pattern) and all(b in IUPAC[c] for c, b in zip(pattern, read))


# ---------------------------------------------------------------------------------------------- generic checks
def deadline_allowed(sc: dict, ans: dict) -> bool:
    """Whether |ans| is a 503 `deadline` of /pattern or /resolve that the scenario's status expectation allows
    implicitly: an answer that could not be written within its time budget is a 503 deadline whatever the budget
    (SPEC-pattern-search §6, §7.6; SPEC-labeled-traversal-core §4.5: never a partial answer), so on a loaded host any
    such request may get it in place of its 200 (the README's `machine-dependent`: time budgets give 200 or 503).
    Allowed where the status expectation asks for a 200 and does not itself name the deadline (an http 503, the code
    `deadline` or a text naming it): a scenario that names it is judged as it states. A deadline before its budget
    fails budget_respected (_budget_503). /resolve only with a request budget (A2-01 of the review of lane S5b):
    SPEC-labeled-traversal-core §4.5 makes its deadline opt-in ("a request without bounds.time_budget_ms runs
    without a deadline ... whatever the server's cap"; capabilities resolve.time_budget.default null), so a 503
    deadline of a /resolve without bounds.time_budget_ms is a host defect: the status fails it."""
    if sc.get("route") not in ("POST /pattern", "POST /resolve") or ans.get("http") != 503 or \
            (ans.get("error") or {}).get("code") != "deadline":
        return False
    if sc.get("route") == "POST /resolve" and _resolve_budget(sc) is None:
        return False
    for exp in sc.get("expect", []):
        if exp.get("check") != "status":
            continue
        args = exp.get("args") or {}
        want = args.get("http", 200)
        wants = want if isinstance(want, list) else [want]
        if 200 not in wants or 503 in wants or "deadline" in json.dumps(args):
            return False
    return True


# the name of before (lane S5, /pattern only)
pattern_deadline_allowed = deadline_allowed

# The front ingress of the public URL (nginx at https://metagraph.ethz.ch:8080, the direct-remote arm; README) cuts a
# request at about FRONT_INGRESS_CUT_MS and answers its own 504, with no host body. The rule (ingress_cut_allowed)
# covers a request whose stated budget B is at least FRONT_INGRESS_CUT_MS - FRONT_INGRESS_CUT_TOLERANCE_MS
# (110,000 ms) and whose 504 came no earlier than that: v6's /pattern budget of 110,000 ms (the owner, 2026-10-16;
# 120,000 ms before, where the cut came before the host's answer could). At B = 110,000 ms the host answers at the
# budget (a 200 limited by it, or its 503 deadline at ~B + 250 ms), about 10 s before the cut, so a 504 at the cut
# says the host had not answered 10 s past B: within the window in which a late 503 still reads inconclusive
# (_budget_503: up to max(1 s, B) past B, here to 220 s), and the 504 hides whether its answer would have come in it
# (inconclusive) or after it (a finding). Nothing the host did is observed: the 503 deadline's case (inconclusive),
# never a host defect and never a pass. The tolerance is the cut's measurement error (~120 s, not stated) and keeps
# v6's budget covered inclusively at B = 110,000 ms; the late window must reach past the cut with it, B + max(1 s, B)
# >= FRONT_INGRESS_CUT_MS + FRONT_INGRESS_CUT_TOLERANCE_MS (test_builder_v6). A request whose budget ends earlier
# (every /resolve and /traverse budget, 60,000 ms, and v5's /pattern ones: the late window ends at the cut) must be
# answered by the host before the cut, so its 504 stays the status failure it was.
FRONT_INGRESS_CUT_MS = 120000.0
FRONT_INGRESS_CUT_TOLERANCE_MS = 10000.0
# the arm behind that ingress
INGRESS_ARM = "direct-remote"
# what time_limited names the ingress's cut (compare._own_time_limitation reads the line's recorded time_limited)
INGRESS_CUT_LIMITATION = "the front ingress's 504 at its cut, past the time budget"


def _stated_budget(sc: dict) -> Optional[float]:
    """The time budget the scenario's request states (ms): /pattern's time_budget_ms, /resolve's
    bounds.time_budget_ms; None when the route or the request names none."""
    route = sc.get("route")
    req = sc.get("request") if isinstance(sc.get("request"), dict) else {}
    if route == "POST /pattern":
        return _number(req.get("time_budget_ms"))
    if route == "POST /resolve":
        return _resolve_budget(sc)
    return None


def ingress_cut_allowed(sc: dict, ans: dict, client: Optional[dict] = None, arm: Optional[str] = None) -> bool:
    """Whether |ans| is the front ingress's 504 in place of the host's answer, which the scenario's status
    expectation allows implicitly: an http 504 without a host body (no error code), on a POST /pattern or
    POST /resolve whose request states a time budget that reaches the ingress's cut (budget >= FRONT_INGRESS_CUT_MS
    - FRONT_INGRESS_CUT_TOLERANCE_MS), delivered at the cut (the client's wall time, when the runner gives it, at
    least that far too), where the status expectation asks for a 200 and does not itself name a 504. The host's
    answer to such a request comes at its budget (SPEC-pattern-search §7.6, SPEC-labeled-traversal-core §4.5) or,
    a 503 deadline, later within the window the checks read inconclusive (_budget_503: max(1 s, the budget) past
    it), and the cut falls in that window, so the 504 shows nothing of the host (as an allowed 503 deadline shows
    nothing: inconclusive). A request whose budget ends before the cut, or a 504 well before it, is not covered, nor
    an arm the ingress does not sit in front of (|arm| given and not direct-remote: direct-local reaches the
    container's port, the MCP's answers are the service's own)."""
    if ans.get("http") != 504 or (ans.get("error") or {}).get("code") is not None:
        return False
    if arm is not None and arm != INGRESS_ARM:
        return False
    budget = _stated_budget(sc)
    if budget is None or budget < FRONT_INGRESS_CUT_MS - FRONT_INGRESS_CUT_TOLERANCE_MS:
        return False
    wall = _number((client or {}).get("wall_ms"))
    if wall is not None and wall < FRONT_INGRESS_CUT_MS - FRONT_INGRESS_CUT_TOLERANCE_MS:
        return False
    for exp in sc.get("expect", []):
        if exp.get("check") != "status":
            continue
        args = exp.get("args") or {}
        want = args.get("http", 200)
        wants = want if isinstance(want, list) else [want]
        if 200 not in wants or 504 in wants or "504" in json.dumps(args):
            return False
    return True


def _resolve_budget(sc: dict) -> Optional[float]:
    """The bounds.time_budget_ms a /resolve scenario's request names (a number), else None."""
    req = sc.get("request") if isinstance(sc.get("request"), dict) else {}
    bounds = req.get("bounds") if isinstance(req.get("bounds"), dict) else {}
    v = bounds.get("time_budget_ms")
    return float(v) if isinstance(v, (int, float)) and not isinstance(v, bool) else None


def _effective_budget(sc: dict, body: Any) -> Optional[float]:
    """The time budget (ms) an answer is held to: its echo when that is a number (/pattern and /resolve
    limits.time_budget_ms, /traverse strategy.bounds.time_budget_ms), else the request's (_stated_budget; /traverse
    strategy.bounds.time_budget_ms); None when neither states one. An echo that is null or no number never hides the
    request's budget (SD-286)."""
    b = body if isinstance(body, dict) else {}
    if sc.get("route") == "POST /traverse":
        def bounds(x: Any) -> dict:
            st = x.get("strategy") if isinstance(x, dict) and isinstance(x.get("strategy"), dict) else {}
            return st.get("bounds") if isinstance(st.get("bounds"), dict) else {}
        echo = _number(bounds(b).get("time_budget_ms"))
        return echo if echo is not None else _number(bounds(sc.get("request")).get("time_budget_ms"))
    lim = b.get("limits") if isinstance(b.get("limits"), dict) else {}
    echo = _number(lim.get("time_budget_ms"))
    return echo if echo is not None else _stated_budget(sc)


def check_status(sc, ans, a, arm, client=None):
    http = ans.get("http")
    want = a.get("http", 200)
    wants = want if isinstance(want, list) else [want]
    if http not in wants:
        err = ans.get("error") or {}
        if deadline_allowed(sc, ans):
            return INCONCLUSIVE, "http 503 deadline in place of %s: the answer could not be written within its " \
                                 "time budget, which %s allows on any host (%s)" % (
                                     want, "SPEC-pattern-search §6/§7.6" if sc.get("route") == "POST /pattern"
                                     else "SPEC-labeled-traversal-core §4.5", (err.get("text") or "")[:200])
        if ingress_cut_allowed(sc, ans, client, arm):
            wall = _number((client or {}).get("wall_ms"))
            return INCONCLUSIVE, "http 504 from the front ingress in place of %s%s: the ingress cuts a request at " \
                                 "~%.0f s, inside the window in which the host's answer to the %.0f ms budget may " \
                                 "still come (its 503 deadline up to max(1 s, the budget) past it reads " \
                                 "inconclusive), so nothing the host did is observed (an allowed deadline's case; " \
                                 "the README's direct-remote)" % (
                                     want, (" at a client wall time of %.0f ms" % wall) if wall is not None else "",
                                     FRONT_INGRESS_CUT_MS / 1000.0, _stated_budget(sc))
        if sc.get("route") == "POST /resolve" and http == 503 and err.get("code") == "deadline" and \
                _resolve_budget(sc) is None:
            return FAIL, "http 503 deadline, expected %s: the request names no bounds.time_budget_ms, and the " \
                         "/resolve deadline is opt-in (SPEC-labeled-traversal-core §4.5: no deadline without the " \
                         "field, whatever the server's cap) (%s)" % (want, (err.get("text") or "")[:200])
        return FAIL, "http %s, expected %s (error %s: %s)" % (http, want, err.get("code"), (err.get("text") or "")[:300])
    err = ans.get("error") or {}
    if "code" in a:
        if err.get("code") != a["code"]:
            return FAIL, "code %r, expected %r (text %r)" % (err.get("code"), a["code"], (err.get("text") or "")[:300])
    text = err.get("text") or ""
    if "text" in a and text != a["text"]:
        return FAIL, "text %r, expected %r" % (text[:400], a["text"])
    if "text_prefix" in a and not text.startswith(a["text_prefix"]):
        return FAIL, "text %r does not start with %r" % (text[:400], a["text_prefix"])
    for frag in a.get("text_contains", []):
        if frag not in text:
            return FAIL, "text %r lacks %r" % (text[:400], frag)
    return PASS, "http %s%s" % (http, (" code %s" % err.get("code")) if err.get("code") else "")


# What the MCP arm's /pattern answer states on its entries in place of the host's top level: the service's
# pattern_search keeps each database's predicate block (normal form, names, known, unknown_labels, vacuous, scope,
# strands) and its timing.load_ms on every entry it merges (by_database[].predicate / .load_ms), and the runner's
# normaliser leaves them there (its body is mode, patterns and index). {the host's top-level path: the entry's key}
MCP_PATTERN_ENTRY_FIELDS = {"predicate": "predicate", "timing.load_ms": "load_ms"}


class _Differs(Exception):
    """The answered entries state a field of MCP_PATTERN_ENTRY_FIELDS differently."""


def _answer_field(sc: dict, body: Any, path: str, arm: str) -> Any:
    """The value of |path| in an answer's body (_get); on the MCP arm, for a POST /pattern path of
    MCP_PATTERN_ENTRY_FIELDS that the body does not state, the value every answered entry states there (one value, the
    host's top-level one: _Differs when two entries differ). Missing when neither states it."""
    try:
        return _get(body, path)
    except Missing:
        if arm != "mcp" or sc.get("route") != "POST /pattern" or not isinstance(body, dict):
            raise
        top = next((k for k in MCP_PATTERN_ENTRY_FIELDS if path == k or path.startswith(k + ".")), None)
        if top is None:
            raise
        sub = MCP_PATTERN_ENTRY_FIELDS[top] + path[len(top):]
        vals = [_get(e, sub) for e in body.get("patterns") or []
                if isinstance(e, dict) and "counts" in e and "error" not in e]
        if not vals:
            raise
        if any(v != vals[0] for v in vals[1:]):
            raise _Differs("%s: the entries state %s" % (path, ", ".join(repr(v) for v in vals[:3])))
        return vals[0]


def check_answer_fields(sc, ans, a, arm):
    """Fields of the answer's body by dotted path (any_status: also on an error body, e.g. a 409's state). when_http:
    the fields of an answer of that status only (another status: nothing to check, pass), e.g. the shape of an
    attempt's duplicate 409 (A-10). On the MCP arm a /pattern field the service states on every entry instead
    (MCP_PATTERN_ENTRY_FIELDS) is read there (_answer_field)."""
    if "when_http" in a and ans.get("http") != a["when_http"]:
        return PASS, "http %s: the fields are checked on an http %s only" % (ans.get("http"), a["when_http"])
    body = _body(ans)
    try:
        for path, want in (a.get("equals") or {}).items():
            got = _answer_field(sc, body, path, arm)
            if got != want:
                return FAIL, "%s = %r, expected %r" % (path, got, want)
        for path in a.get("present", []):
            _answer_field(sc, body, path, arm)
        for path in a.get("absent", []):
            try:
                _answer_field(sc, body, path, arm)
                return FAIL, "%s present, expected absent" % path
            except Missing:
                pass
    except _Differs as d:
        return FAIL, "%s (one database, one value)" % d
    return PASS, "fields as expected"


def check_entry_fields(sc, ans, a, arm):
    """The entry's fields: `equals` (a float within `tolerance`), `present`, `absent`. Every field the answer carries
    is compared before a missing one is reported (Missing: fail on a direct arm, not exposed on the MCP's), so a
    missing field never hides a wrong one beside it (the MCP's entries carry `scope` but not the host's absence_scope:
    a mark that sets the missing field aside never sets aside a failure; VERIFY of round 2, L3)."""
    e = _entry(ans, a.get("pattern", 0))
    missing = []
    for path, want in (a.get("equals") or {}).items():
        try:
            got = _get(e, path)
        except Missing:
            missing.append(path)
            continue
        if isinstance(want, float) and isinstance(got, (int, float)):
            # information bits and the like are sums of exact binary logarithms: equal up to rounding only
            if abs(got - want) > a.get("tolerance", 1e-9):
                return FAIL, "%s = %r, expected %r" % (path, got, want)
        elif got != want:
            return FAIL, "%s = %r, expected %r" % (path, got, want)
    for path in a.get("present", []):
        try:
            _get(e, path)
        except Missing:
            missing.append(path)
    for path in a.get("absent", []):
        try:
            _get(e, path)
            return FAIL, "%s present, expected absent" % path
        except Missing:
            pass
    if missing:
        raise Missing(missing[0])
    return PASS, "entry fields as expected"


def _time_limited(ans: dict) -> str:
    """What makes an answer time-limited (SPEC-labeled-traversal-core §6.8: determinism holds provided no time
    budget tripped; SPEC-pattern-search §7.9: determinism "time_limited", and "a 503 deadline depends on the machine
    too"); '' when nothing did. A deadline answer first: the host's 503 deadline (/pattern, /resolve) and its relays
    by the service (a failed search whose error starts pattern_deadline, traversal_resolve_deadline)."""
    err = ans.get("error") if isinstance(ans.get("error"), dict) else {}
    http, code = ans.get("http"), err.get("code")
    if http == 503 and code == "deadline":
        return "a 503 deadline: nothing written within the time budget"
    if http == 502 and code == "search_failed" and str(err.get("text") or "").startswith("pattern_deadline"):
        return "the host's /pattern 503 deadline relayed by the service"
    if http == 503 and code == "traversal_resolve_deadline":
        return "the host's /resolve 503 deadline relayed by the service"
    body = ans.get("body")
    if not isinstance(body, dict):
        return ""
    for i, e in enumerate(body.get("patterns") or []):
        if isinstance(e, dict) and (e.get("determinism") == "time_limited" or (e.get("stop") or {}).get("reason") == "time"
                                    or (e.get("cut") or {}).get("reason") == "time"):
            return "pattern %d time-limited" % i
    if (body.get("stop") or {}).get("reason") == "time":
        return "resolve stopped by time"
    for i, r in enumerate(body.get("results") or []):
        if not isinstance(r, dict):
            continue
        if (r.get("resource_stop") or {}).get("resource") == "time":
            return "seed %d stopped by its time budget" % i
        # the attempt's duration bound: the seeds' time budgets plus an allowance (SPEC §6.8)
        if (r.get("resource_stop") or {}).get("resource") == "attempt_deadline":
            return "seed %d stopped by its attempt's duration bound" % i
        if "arms" not in r and _AT_BOUND_RE.fullmatch(_seed_error_text(r) or ""):
            return "seed %d: the attempt reached its duration bound (AttemptAtBound)" % i
        # a derivation the time budget ended (failed, or partial: traverse.cpp failed_seed_to_json, the D3 partial)
        for lim in r.get("limitations") if isinstance(r.get("limitations"), list) else []:
            if isinstance(lim, dict) and lim.get("knob") == _TIME_KNOB:
                return "seed %d: its %s limitation names the time budget" % (i, lim.get("kind"))
        for side, armo in (r.get("arms") or {}).items():
            if (armo.get("cap_trigger") or {}).get("reason") == "time_budget":
                return "seed %d %s arm stopped by time_budget" % (i, side)
            # a time cap that cut walks after another trigger (traverse.cpp arm_limitations: one walk_domain
            # limitation per cap that ended walks)
            for lim in armo.get("limitations") if isinstance(armo.get("limitations"), list) else []:
                if isinstance(lim, dict) and lim.get("knob") == _TIME_KNOB:
                    return "seed %d %s arm: its %s limitation names the time budget" % (i, side, lim.get("kind"))
    return ""


# the knob a /traverse limitation names when a time budget ended the work it states (a walk cap or a derivation)
_TIME_KNOB = "bounds.time_budget_ms"


def check_deterministic(sc, ans, a, arm, previous=None, client=None):
    """The answer's digest equals the earlier repeats' (`previous`: their digests, or {"digest", "time_limited"}
    objects; a bare digest counts as not time-limited). Only earlier repeats no time budget limited are compared (a
    time-limited answer, a 503 deadline included, need not repeat: SPEC-pattern-search §7.9,
    SPEC-labeled-traversal-core §6.8): with none left the check is inconclusive, and an answer that is itself
    time-limited is inconclusive unless its stated time limitation cannot have happened (fail;
    _time_limitation_problem, a /resolve stop judged against the client's wall, |client|)."""
    d = answer_digest(ans, sc.get("route"))
    if not previous:
        return INCONCLUSIVE, "first repeat: nothing to compare yet (digest %s)" % d[:16]
    prev = [{"digest": p.get("digest"), "time_limited": p.get("time_limited") or ""} if isinstance(p, dict)
            else {"digest": p, "time_limited": ""} for p in previous]
    own = time_limited(ans, sc, client, arm)
    clean = [p for p in prev if not p["time_limited"]]
    diff = [p for p in clean if p["digest"] != d]
    if own:
        if diff:
            bad = _time_limitation_problem(sc, ans, arm, client)
            if bad:
                return FAIL, "digest %s differs from %d earlier repeat(s); its time limitation cannot have happened: " \
                             "%s" % (d[:16], len(diff), bad)
        return INCONCLUSIVE, "a time budget limited this repeat (%s): it need not repeat (digest %s%s)" % (
            own, d[:16], ", differs from %d earlier repeat(s) without one" % len(diff) if diff else "")
    if not clean:
        return INCONCLUSIVE, "no earlier repeat without a time limitation (%s): nothing to compare (digest %s)" % (
            prev[0]["time_limited"], d[:16])
    if diff:
        return FAIL, "digest %s differs from %d earlier repeat(s) without a time limitation (%s)" % (
            d[:16], len(diff), str(diff[0]["digest"])[:16])
    skipped = len(prev) - len(clean)
    return PASS, "digest %s equal to %d earlier repeat(s)%s" % (
        d[:16], len(clean), " (%d time-limited one(s) not compared)" % skipped if skipped else "")


def check_content_encoding(sc, ans, a, arm):
    h = ans.get("headers")
    if h is None:
        raise Missing("headers")
    got = h.get("content-encoding")
    if got != a.get("value"):
        return FAIL, "content-encoding %r, expected %r" % (got, a.get("value"))
    return PASS, "content-encoding %r" % got


def check_capabilities_doc(sc, ans, a, arm):
    body = _body(ans)
    for path in a.get("required", []):
        _get(body, path)
    for path, want in (a.get("equals") or {}).items():
        got = _get(body, path)
        if got != want:
            return FAIL, "%s = %r, expected %r" % (path, got, want)
    if "max_bytes" in a:
        n = len(json.dumps(body, separators=(",", ":")).encode())
        if n > a["max_bytes"]:
            return FAIL, "document of %d bytes, above %d" % (n, a["max_bytes"])
    return PASS, "capabilities document as expected"


# ---------------------------------------------------------------------------------------------- stated reasons
# A check that cannot test its property because the answer states a hiding (a withheld list, a cut), a stop or a
# time limitation accepts the statement only once the answer shows it can have happened: a stated reason that
# cannot have happened FAILS the check that met it (and the route's invariants), it never makes it inconclusive.
# The reasons and what must hold for each are listed in README.md, "Stated reasons"; in short, /pattern
# (SPEC-pattern-search §7.5, §7.6): a threshold reason needs its count above its threshold, a budget reason its
# budget reached (max_steps: the steps of all patterns at max_steps), a time reason a work time that had passed,
# a stop_at_threshold stop the knob set, a label or predicate reason its stop or its refused rows; /resolve
# (SPEC-labeled-traversal-core §4.5): a time stop at the budget less the reserve; /traverse (§6.7, §7.0): a cap
# trigger the walk_domain limitation of its knob with the observed demand past the limit, the arm's counters at
# the cap, a time stop at or after the budget.

# SPEC-pattern-search §7.6: E = 1.25 x (B / (b x 1000) + B / (c x 1000)) ms for B bytes of results built, at the
# delivery rates b and c in MB/s (--pattern-delivery-build-mbps 10, --pattern-delivery-compress-mbps 50 by default).
_DELIVERY_MBPS = (10.0, 50.0)
# The finalisation reserves (ms): /resolve's is the host's constant (SPEC-labeled-traversal-core §4.5; capabilities
# resolve.time_budget.finalize_reserve_ms, finalize_reserve_configurable false), so a reserve the answer states is
# checked against it, never believed (suite_common's rule: never from the failing answer); /pattern's is the host's
# configured one (SPEC-pattern-search §10.2, --pattern-finalize-ms, default 250), which only the answer states here
_RESOLVE_RESERVE_MS = 250.0
_PATTERN_RESERVE_MS = 250.0

# Bytes counted for each result an answer withholds (built, then not handed over): generous, so that E is an over-
# estimate and the work time a lower bound (an occurrence list makes labelled results far larger).
_RESULT_BYTES = {"plain": 2048, "labelled": 65536}


def _limit(sc: dict, body: dict, name: str) -> Any:
    """The effective value of a /pattern knob: the answer's `limits`, else the request's."""
    lim = body.get("limits") if isinstance(body, dict) else None
    if isinstance(lim, dict) and name in lim and lim[name] is not None:
        return lim[name]
    req = sc.get("request") if isinstance(sc.get("request"), dict) else {}
    return req.get(name)


def _pattern_work_time(sc: dict, body: dict) -> Optional[float]:
    """A lower bound of the work time of SPEC-pattern-search §7.6, time_budget_ms - finalize_reserve_ms - E: E is
    taken at twice its formula at the default (slowest) delivery rates, with every result the answer withholds
    counted at a generous size (_RESULT_BYTES), so an elapsed time below this bound is a time stop or a time
    limitation that came before the work time had passed. None without a budget."""
    budget = _effective_budget(sc, body)
    if budget is None:
        return None
    reserve = _number((body.get("limits") or {}).get("finalize_reserve_ms"))
    reserve = _PATTERN_RESERVE_MS if reserve is None else reserve
    req = sc.get("request") if isinstance(sc.get("request"), dict) else {}
    labelled = (req.get("output") or {}).get("labels") in ("all", "predicate_only") or req.get("predicate") is not None
    per = _RESULT_BYTES["labelled" if labelled else "plain"]
    nbytes = 0
    for e in body.get("patterns") or []:
        if not isinstance(e, dict) or "counts" not in e:
            continue
        nbytes += len(json.dumps(e.get("results") or [], separators=(",", ":")))
        if e.get("withheld"):
            long_ = "contexts" not in (e.get("counts") or {})
            cap = _number(_limit(sc, body, "max_paths" if long_ else "max_contexts")) or 0
            c = (e.get("counts") or {}).get(_primary_count_name(e, req)) or {}
            n = c.get("value") if c.get("relation") in ("exact", "at_least") and isinstance(c.get("value"), int) else cap
            nbytes += int(min(n, cap)) * per
    b, c = _DELIVERY_MBPS
    e_ms = 2 * 1.25 * (nbytes / (b * 1000.0) + nbytes / (c * 1000.0))
    return budget - reserve - e_ms


# The request's time outside its patterns' own timing (the gaps between the patterns' searches): an answer judged by
# its entries' timing (_pattern_elapsed) may have run this much longer than they state. At most 4.6 ms over the 974
# direct /pattern 200s of the v5 baseline run (top-level elapsed less the entries' sum, never below 0).
_ENTRY_TIMING_SLACK_MS = 25.0
# An entry's timing values that add up to its part of the request's time (SPEC-pattern-search §8.4): elapsed_ms, that
# pattern's search, and the reads and passes it states beside the search (label_discovery_ms, placement_ms,
# selection_ms). The others are parts of its search (extension_ms, label_intersection_ms, verification_ms of a long
# pattern's search; support_ms of extension_ms) and are not added: summed, they ran up to 1,047 ms past the top-level
# figure over the v5 baseline's direct answers (L-long-0003), which would admit a stop that early
_ENTRY_TIMING_OWN = ("elapsed_ms", "label_discovery_ms", "placement_ms", "selection_ms")


def _pattern_elapsed(body: dict) -> Tuple[Optional[float], str]:
    """(the elapsed time a /pattern answer states for its search, in ms, where it comes from): the top-level
    timing.elapsed_ms (from the deadline's start to the end of the search, SPEC-pattern-search §8.4); else, on an
    answer without it (the MCP arm's: the service keeps each pattern's own timing and not the request's), the sum
    over the answered entries of the timing values _ENTRY_TIMING_OWN names (elapsed_ms is that pattern's search; the
    reads beside it, label_discovery_ms, placement_ms, selection_ms, are its own values; the patterns are searched
    one after the other within the request's deadline), plus _ENTRY_TIMING_SLACK_MS; (None, '') when neither states
    a time."""
    el = _number((body.get("timing") or {}).get("elapsed_ms"))
    if el is not None:
        return el, "timing.elapsed_ms"
    total, seen = 0.0, 0
    for e in body.get("patterns") or []:
        if not isinstance(e, dict) or "counts" not in e or "error" in e:
            continue
        t = e.get("timing") if isinstance(e.get("timing"), dict) else {}
        if _number(t.get("elapsed_ms")) is None:
            return None, ""
        seen += 1
        total += sum(_number(t.get(k)) or 0.0 for k in _ENTRY_TIMING_OWN)
    if not seen:
        return None, ""
    return total + _ENTRY_TIMING_SLACK_MS, "the entries' timing (%d pattern(s), +%.0f ms)" % (
        seen, _ENTRY_TIMING_SLACK_MS)


def _pattern_time_problem(sc: dict, body: dict, what: str) -> Optional[str]:
    """Why a /pattern time stop or time limitation (|what|) cannot have happened, or None: the work stops at the work
    time and the answer's elapsed time runs to the end of the search (SPEC-pattern-search §7.6), the top-level
    timing.elapsed_ms, else what its entries state (_pattern_elapsed: the MCP arm's answer has no top-level timing).
    An answer that states no time at all cannot show it: Undecidable (the check that met the statement reads
    inconclusive, never pass)."""
    w = _pattern_work_time(sc, body)
    if w is None:
        return None
    el, source = _pattern_elapsed(body)
    if el is None:
        raise Undecidable("%s, and the answer states no timing (neither timing.elapsed_ms nor an entry's): it cannot "
                          "be held to the work time (at least %.0f ms of the %s ms budget: SPEC-pattern-search §7.6)" % (
                              what, w, _limit(sc, body, "time_budget_ms")))
    if el >= w:
        return None
    return "%s after %.0f ms (%s), before the work time (at least %.0f ms of the %s ms budget: SPEC-pattern-search " \
           "§7.6)" % (what, el, source, w, _limit(sc, body, "time_budget_ms"))


def _count_above(c: Any, cap: Optional[float]) -> Optional[bool]:
    """Whether a count is shown above |cap| (exact > cap, or bounds whose upper bound > cap); None when the count
    cannot show it either way (at_least, unknown) or the cap is unknown."""
    if not isinstance(c, dict) or cap is None:
        return None
    if c.get("relation") == "exact" and isinstance(c.get("value"), (int, float)):
        return c["value"] > cap
    if c.get("relation") == "bounds" and isinstance(c.get("upper"), (int, float)):
        return c["upper"] > cap
    if c.get("relation") == "at_least" and isinstance(c.get("value"), (int, float)) and c["value"] > cap:
        return True
    return None


def _count_text(c: Any) -> str:
    if not isinstance(c, dict):
        return "no count"
    if c.get("relation") == "bounds":
        return "bounds %s..%s" % (c.get("lower"), c.get("upper"))
    return "%s %s" % (c.get("relation"), c.get("value"))


def _pattern_reason_problems(sc: dict, body: dict, i: int, e: dict) -> List[str]:
    """What is wrong with entry i's stated withheld reason, cut reason, stop and time limitation: each must be one
    SPEC-pattern-search §7.5 / §7.6 lists, possible for the request (its mode, whether it reads annotation, has a
    memory account, a predicate, a path search) and must have happened by what the answer shows (README, "Stated
    reasons" and "Lane S6g1"). Nothing is judged where the answer does not show the quantity (a missing cap, work or
    timing)."""
    probs: List[str] = []
    if not isinstance(e, dict) or "counts" not in e or "error" in e:
        return probs
    req = sc.get("request") if isinstance(sc.get("request"), dict) else {}
    counts = e.get("counts") or {}

    def cnt(name: str) -> Optional[dict]:
        c = counts.get(name)
        return c if isinstance(c, dict) else None
    long_ = "contexts" not in counts
    main = cnt("contexts") if not long_ else (
        cnt("supported_paths") if req.get("long_search") == "supported_paths" and cnt("supported_paths") else cnt("paths"))
    main_cap = _number(_limit(sc, body, "max_paths" if long_ else "max_contexts"))
    raw = cnt("contexts") or cnt("supported_paths") or cnt("paths")
    stop = e.get("stop") if isinstance(e.get("stop"), dict) else {}
    sr = stop.get("reason")
    threshold = _limit(sc, body, "stop_at_threshold") is True
    predicate = req.get("predicate") is not None
    labels = (req.get("output") or {}).get("labels")
    time_limited = sr == "time" or e.get("determinism") == "time_limited"
    max_steps = _number(_limit(sc, body, "max_steps"))
    steps = [((x.get("work") or {}).get("steps")) for x in body.get("patterns") or []
             if isinstance(x, dict) and "counts" in x]
    steps_sum = sum(steps) if steps and all(isinstance(v, int) for v in steps) else None
    steps_short = max_steps is not None and steps_sum is not None and steps_sum < max_steps
    where = "p%d" % i
    # what the request reads and holds: annotation rows are read only with labels all in a retrieval mode (§14.1:
    # mode count reads none whatever the projection), a supported-path search (§20.3) or a predicate (§19.6, in
    # every mode); the label-free path (mode count, or labels none) has a memory account only for a path search, and
    # labels all in a retrieval mode and a predicate add one (§7.6 "Memory", §14.4, §19.9)
    mode = e.get("mode", req.get("mode", "all_or_count"))
    ls = req.get("long_search")
    path_search = long_ and ls in ("paths", "supported_paths")
    labels_read = labels == "all" and mode != "count"
    reads_annotation = labels_read or (long_ and ls == "supported_paths") or predicate
    memory_account = labels_read or predicate or path_search
    nothing_read = "a request that reads no annotation (labels %r in mode %s, no predicate, long_search %r)" % (
        labels, mode, ls)
    no_account = "a request without a memory account (labels %r in mode %s, no predicate, no path search: " \
                 "§7.6)" % (labels, mode)
    if e.get("rows_refused") and not reads_annotation:
        probs.append("%s: rows_refused on %s" % (where, nothing_read))
    # the request's selection work as the answer shows it, through this pattern (§19.9: one budget per request,
    # the pass stops once the units reached max_predicate_work; a later pattern finds it spent, sticky: §19.8,
    # §20.9); None where the answer does not show it
    pmax = _number(_limit(sc, body, "max_predicate_work"))
    punits = [_number((x.get("work") or {}).get("predicate_units")) for x in (body.get("patterns") or [])[:i + 1]
              if isinstance(x, dict) and "counts" in x and "error" not in x]
    pwork_spent = None if pmax is None or not punits or any(u is None for u in punits) else sum(punits) >= pmax
    own_units = _number((e.get("work") or {}).get("predicate_units"))

    w = e.get("withheld")
    wr = w.get("reason") if isinstance(w, dict) else (w if isinstance(w, str) else None)
    if wr not in (None, "paths_later_increment", "anchors_above_threshold") and mode != "all_or_count":
        # only all_or_count withholds, but for these two of patterns longer than k (§7.5, §12.1)
        probs.append("%s: withheld %r in mode %s" % (where, wr, mode))
    if wr == "count_above_threshold":
        # exact (or, without the mask, bounds) and above the threshold: a lower bound is threshold_crossed's
        constant_true = (e.get("selection") or {}).get("pass") == "constant" and \
            (body.get("predicate") or {}).get("normal_form") is True
        if predicate and not long_ and not constant_true:
            # a normal form that is the constant `true` runs the unfiltered request with the engine's own thresholds
            # (§19.6: `true` is the unfiltered request), whose count_above_threshold stands; `false` selects nothing
            # and withholds nothing
            probs.append("%s: withheld count_above_threshold with a predicate on L <= k (§19.8: its raw threshold "
                         "is predicate_above_threshold)" % where)
        elif main_cap is not None and not ((main or {}).get("relation") in ("exact", "bounds") and
                                           _count_above(main, main_cap)):
            probs.append("%s: withheld count_above_threshold with %s, not above the threshold %s" % (
                where, _count_text(main), int(main_cap)))
    elif wr == "threshold_crossed":
        if not threshold:
            probs.append("%s: withheld threshold_crossed without stop_at_threshold" % where)
        elif sr not in _THRESHOLD_STOPS:
            probs.append("%s: withheld threshold_crossed without a stop_at_threshold stop (stop %r)" % (where, sr))
    elif wr == "discovery_budget":
        if sr != "max_steps":
            probs.append("%s: withheld discovery_budget without a max_steps stop (stop %r)" % (where, sr))
        elif steps_short:
            probs.append("%s: withheld discovery_budget after %d steps of max_steps %d" % (where, steps_sum, max_steps))
        elif (main or {}).get("relation") == "exact":
            probs.append("%s: withheld discovery_budget with %s: a count the stop interrupted is not exact (§7.4)" % (
                where, _count_text(main)))
    elif wr == "deadline":
        if not time_limited:
            probs.append("%s: withheld deadline without a time stop or time limitation (stop %r)" % (where, sr))
    elif wr == "paths_later_increment":
        if not long_ or ls in ("paths", "supported_paths"):
            probs.append("%s: withheld paths_later_increment for a pattern of L <= k or with long_search %r" % (
                where, ls))
        elif (cnt("anchors") or {}).get("relation") == "exact" and (cnt("anchors") or {}).get("value") == 0:
            probs.append("%s: withheld paths_later_increment with the anchors exact 0 (§7.7: that empty answer is "
                         "complete)" % where)
    elif wr == "anchors_above_threshold":
        cap = _number(_limit(sc, body, "max_anchors"))
        if not path_search:
            probs.append("%s: withheld anchors_above_threshold without a path search (long_search %r: §12.1)" % (
                where, ls))
        elif cap is not None and not _count_above(cnt("anchors"), cap):
            probs.append("%s: withheld anchors_above_threshold with anchors %s, not above max_anchors %s" % (
                where, _count_text(cnt("anchors")), int(cap)))
    elif wr == "annotation_budget":
        if not reads_annotation:
            probs.append("%s: withheld annotation_budget on %s" % (where, nothing_read))
        elif sr not in ("max_annotation_work", "max_memory") and not e.get("rows_refused"):
            probs.append("%s: withheld annotation_budget without an annotation stop or a refused row" % where)
    elif wr == "anchor_labels_truncated":
        if labels != "all":
            probs.append("%s: withheld anchor_labels_truncated without labels all (labels %r: §14.6)" % (where, labels))
        elif not e.get("anchors_truncated"):
            probs.append("%s: withheld anchor_labels_truncated without anchors_truncated" % where)
    elif wr == "output_budget":
        if not memory_account:
            probs.append("%s: withheld output_budget on %s" % (where, no_account))
        elif sr != "max_memory" and e.get("labels_status") != "output_budget":
            probs.append("%s: withheld output_budget without a memory stop" % where)
    elif wr == "predicate_above_threshold":
        cap = _number(_limit(sc, body, "max_predicate_contexts"))
        if not predicate:
            probs.append("%s: withheld predicate_above_threshold without a predicate" % where)
        elif cap is not None and not any(_count_above(cnt(n), cap) for n in ("contexts", "paths", "supported_paths")):
            probs.append("%s: withheld predicate_above_threshold with %s, not above max_predicate_contexts %s" % (
                where, _count_text(raw), int(cap)))
    elif wr == "selected_above_threshold":
        cap = _number(_limit(sc, body, "max_paths" if long_ else "max_contexts"))
        if not predicate:
            probs.append("%s: withheld selected_above_threshold without a predicate" % where)
        elif cap is not None and not _count_above(cnt("selected"), cap):
            probs.append("%s: withheld selected_above_threshold with selected %s, not above %s" % (
                where, _count_text(cnt("selected")), int(cap)))
    elif wr == "predicate_budget":
        if not predicate:
            probs.append("%s: withheld predicate_budget without a predicate" % where)
        elif sr not in ("max_predicate_work", "max_memory") and not e.get("rows_refused") and not pwork_spent:
            # the selection work spent by an earlier pattern (sticky, §19.8) withholds whatever stopped this
            # pattern's discovery first (beside an engine stop of its own: the §19.8 corner pending the engine
            # owner's decision, see the cut max_predicate_work rule)
            probs.append("%s: withheld predicate_budget without a predicate stop, a refused row or the request's "
                         "max_predicate_work spent" % where)
    elif wr is not None:
        probs.append("%s: withheld %r, a reason SPEC-pattern-search §7.5 does not list" % (where, wr))

    cu = e.get("cut")
    cr = cu.get("reason") if isinstance(cu, dict) else (cu if isinstance(cu, str) else None)
    if cr == "max_contexts":
        cap = _number(_limit(sc, body, "max_contexts"))
        c = cnt("selected") if predicate and cnt("selected") else main
        if cap is not None and sr != "max_contexts" and _count_above(c, cap) is False:
            probs.append("%s: cut max_contexts with %s, within max_contexts %s" % (where, _count_text(c), int(cap)))
    elif cr == "max_steps":
        if sr != "max_steps":
            probs.append("%s: cut max_steps without a max_steps stop (stop %r)" % (where, sr))
        elif steps_short:
            probs.append("%s: cut max_steps after %d steps of max_steps %d" % (where, steps_sum, max_steps))
    elif cr == "time":
        if not time_limited:
            probs.append("%s: cut time without a time stop or time limitation (stop %r)" % (where, sr))
    elif cr == "max_memory":
        if not memory_account:
            probs.append("%s: cut max_memory on %s" % (where, no_account))
        elif sr is None and not e.get("rows_refused") and e.get("labels_status") != "output_budget":
            # a memory cut of the list (the descriptors, the result objects, a label-free path's object) replaces
            # any other cut and states stop {output, max_memory} unless an earlier stop is stated (§7.5, §12.1,
            # §14.4, §19.8, §20.5: first stop wins), so only a cut without any stop or refused row is impossible.
            # Lenient beside an earlier stop: the cut states no phase of its own (no cut.phase), so a memory cut
            # that never happened cannot be told from one after that stop
            probs.append("%s: cut max_memory without a memory stop, an earlier stop or a refused row" % where)
    elif cr == "max_paths":
        cap = _number(_limit(sc, body, "max_paths"))
        c = cnt("selected") if predicate and cnt("selected") else main
        if not path_search:
            probs.append("%s: cut max_paths without a path search (long_search %r)" % (where, ls))
        elif cap is not None and sr != "max_paths" and _count_above(c, cap) is False:
            probs.append("%s: cut max_paths with %s, within max_paths %s" % (where, _count_text(c), int(cap)))
    elif cr == "max_annotation_work":
        if not (long_ and ls == "supported_paths"):
            # labels all states an annotation stop without a cut (§7.5); only the supported-path search cuts so
            probs.append("%s: cut max_annotation_work without a supported-path search (long_search %r: §20.5)" % (
                where, ls))
        elif sr != "max_annotation_work":
            probs.append("%s: cut max_annotation_work without that stop (stop %r)" % (where, sr))
    elif cr == "max_predicate_contexts":
        cap = _number(_limit(sc, body, "max_predicate_contexts"))
        if not predicate:
            probs.append("%s: cut max_predicate_contexts without a predicate" % where)
        elif cap is not None and sr != "max_predicate_contexts" and _count_above(raw, cap) is False:
            probs.append("%s: cut max_predicate_contexts with %s, within %s" % (where, _count_text(raw), int(cap)))
    elif cr == "max_predicate_work":
        # §19.8: the cut is the first that applies in the order engine stop (max_steps, time), the pass's stop, the
        # raw release's cut, the selected list's cut, while stop names the first stop in time: a stop_at_threshold
        # stop of the raw discovery (max_predicate_contexts) keeps the stop and lets the pass's cut through, and a
        # selection work spent by an earlier pattern cuts whatever stopped this pattern's discovery (sticky, §19.8,
        # §20.9: pattern.cpp put_selection, put_supported; pattern_selection.cpp select)
        if not predicate:
            probs.append("%s: cut max_predicate_work without a predicate" % where)
        elif sr is None:
            probs.append("%s: cut max_predicate_work without that stop or an earlier one" % where)
        elif sr != "max_predicate_work" and pwork_spent is False:
            probs.append("%s: cut max_predicate_work with the request's selection work %d below max_predicate_work "
                         "%d (stop %r)" % (where, sum(punits), int(pmax), sr))
        elif sr in ("max_steps", "time") and stop.get("phase") in ("discovery", "mask_scan", "extraction") and \
                not long_ and (own_units or 0) > 0:
            # this pattern's own pass ran after an engine stop: the engine's cut comes first (§19.8)
            probs.append("%s: cut max_predicate_work after the engine's stop %s/%s, whose cut comes first (§19.8)" % (
                where, stop.get("phase"), sr))
        # Pending the engine owner's decision (§19.8 corner): a selection budget an earlier pattern spent, beside this
        # pattern's own engine stop (max_steps, time; no pass of its own, own_units 0). §19.8's order puts the
        # engine's cut first (cut max_steps / time, the rule of those cuts above); its sticky sentence, and the
        # route (a fresh extraction without the engine's cut: pattern_selection.cpp select's `before` branch), give
        # cut max_predicate_work in partial, withheld predicate_budget in all_or_count. Both are accepted until §19.8
        # or the route states one.
    elif cr == "max_anchors":
        if sr != "max_anchors":
            probs.append("%s: cut max_anchors without that stop (stop %r)" % (where, sr))
    elif cr is not None:
        probs.append("%s: cut %r, a reason SPEC-pattern-search §7.5 does not list" % (where, cr))

    if sr == "max_steps":
        if steps_short:
            probs.append("%s: a max_steps stop after %d steps of max_steps %d (the steps of all patterns reach it: "
                         "SPEC-pattern-search §7.6)" % (where, steps_sum, max_steps))
    elif sr in _THRESHOLD_STOPS:
        stopped = _threshold_stop_count(sr, stop.get("phase"), cnt, main, predicate, long_)
        if not threshold:
            probs.append("%s: a %s stop (stop_at_threshold) without stop_at_threshold" % (where, sr))
        elif sr == "max_predicate_contexts" and not predicate:
            probs.append("%s: a max_predicate_contexts stop without a predicate" % where)
        elif sr == "max_paths" and not path_search:
            probs.append("%s: a max_paths stop without a path search (long_search %r)" % (where, ls))
        elif sr == "max_anchors" and not long_:
            probs.append("%s: a max_anchors stop of a pattern of L <= k" % where)
        elif sr == "max_contexts" and long_ and not predicate:
            probs.append("%s: a max_contexts stop of a pattern longer than k" % where)
        elif (stopped or {}).get("relation") == "exact":
            # §7.4: exact means no step, time or threshold stop touched the count
            probs.append("%s: a %s stop (stop_at_threshold) with %s: a count the stop interrupted is not exact "
                         "(§7.4)" % (where, sr, _count_text(stopped)))
        else:
            m = _threshold_value_problem(sc, body, sr, stop.get("phase"), stopped, predicate, long_)
            if m:
                probs.append("%s: %s" % (where, m))
    elif sr == "max_memory" and not memory_account:
        probs.append("%s: a max_memory stop of %s" % (where, no_account))
    elif sr == "max_annotation_work":
        if not reads_annotation:
            probs.append("%s: a max_annotation_work stop of %s" % (where, nothing_read))
    elif sr == "max_predicate_work":
        if not predicate:
            probs.append("%s: a max_predicate_work stop without a predicate" % where)
        elif pwork_spent is False:
            probs.append("%s: a max_predicate_work stop with the request's selection work %d below %d (§19.9: the "
                         "pass stops once the units reached it)" % (where, sum(punits), int(pmax)))
    elif sr not in (None, "time", "max_memory"):
        probs.append("%s: stop reason %r, one SPEC-pattern-search §7.6 does not list" % (where, sr))
    # where the stop came, by its reason and the request: SPEC-pattern-search §8.7's phases and reasons
    m = _stop_phase_problem(stop.get("phase"), sr, mode=mode, labels=labels, predicate=predicate,
                            motif=req.get("predicate_scope") == "motif", long_=long_, ls=ls,
                            selection_pass=(e.get("selection") or {}).get("pass")) if stop else None
    if m:
        probs.append("%s: %s" % (where, m))
    if time_limited:
        try:
            m = _pattern_time_problem(sc, body, "%s: %s" % (where, "a time stop" if sr == "time" else "a time limitation"))
        except Undecidable:
            if probs:
                return probs   # a stated reason that cannot have happened decides first
            raise
        if m:
            probs.append(m)
    return probs


# the stops of stop_at_threshold (SPEC-pattern-search §7.6, §12.1, §19.8, §20.5, §20.9)
_THRESHOLD_STOPS = ("max_contexts", "max_anchors", "max_paths", "max_predicate_contexts")

# the phases of a stop and the reasons each states (SPEC-pattern-search §7.6, §8.7 `stop`, §12.1, §14.4, §19.8,
# §20.5, §20.9, §25.5; the engine's StopPhase in graph/alignment/pattern_search.hpp, the route's set_stop calls in
# src/cli/pattern*.cpp)
_STOP_PHASE_REASONS = {
    "discovery": {"max_steps", "time", "max_contexts", "max_anchors", "max_predicate_contexts"},
    "mask_scan": {"max_steps", "time"},
    "extraction": {"time"},
    "extension": {"max_steps", "time", "max_paths", "max_memory", "max_annotation_work", "max_predicate_contexts"},
    "label_discovery": {"max_annotation_work", "max_memory", "time"},
    "placement": {"max_annotation_work", "max_memory", "time"},
    "output": {"max_memory", "time"},
    "selection": {"max_predicate_work", "max_memory", "time", "max_contexts", "max_annotation_work"},
}


def _stop_phase_problem(phase: Any, reason: Any, mode: str, labels: Any, predicate: bool, motif: bool, long_: bool,
                        ls: Any, selection_pass: Any) -> Optional[str]:
    """Why a stop {phase, reason} cannot have come in this request, or None: the phase must be one §8.7 lists, state
    that reason, and be a phase the request runs. extraction: a release (a retrieval mode, or a predicate's raw
    release, which count runs too: §19.6); extension: a path search of a pattern longer than k (§12.1, §20.5);
    label_discovery and placement: a projection's reads in a retrieval mode (labels all, predicate_only: §14.1,
    §19.2); output: the results a retrieval mode builds beyond the engine's release (labels, a predicate's selected
    contexts, a path search's list), or a motif's evaluation (§25.5); selection: a predicate (§19.8). By reason:
    max_contexts in discovery only for a pattern of L <= k without a predicate pass (a constant normal form runs the
    engine's threshold; a pass, stated by selection.pass, renames it max_predicate_contexts: pattern.cpp
    put_selection), in selection only for a predicate on L <= k; max_predicate_contexts in discovery only on L <= k,
    in extension only for supported paths (§20.9); max_annotation_work in extension and selection only for supported
    paths (§20.5, §20.9). A stop without a phase is left to the schema."""
    if phase is None:
        return None
    allowed = _STOP_PHASE_REASONS.get(phase)
    if allowed is None:
        return "stop phase %r, one SPEC-pattern-search §8.7 does not list" % (phase,)
    if reason is not None and reason not in allowed:
        return "a %s stop in the phase %s, which states only %s (SPEC-pattern-search §8.7)" % (
            reason, phase, "/".join(sorted(allowed)))
    retrieval = mode != "count"
    path_search = long_ and ls in ("paths", "supported_paths")
    if phase == "extraction" and not (retrieval or predicate):
        return "an extraction stop of mode count without a predicate (nothing released: §7.6)"
    if phase == "extension" and not path_search:
        return "an extension stop without a path search (%s, long_search %r: §12.1)" % (
            "a pattern longer than k" if long_ else "a pattern of L <= k", ls)
    if phase in ("label_discovery", "placement") and not (retrieval and labels in ("all", "predicate_only")):
        return "a %s stop of a request that reads no labels (labels %r in mode %s: §14.1, §19.2)" % (phase, labels, mode)
    if phase == "output" and not ((retrieval and (labels in ("all", "predicate_only") or predicate or path_search))
                                  or (predicate and motif)):
        return "an output stop of a request that builds no output beyond the engine's release (labels %r in mode " \
               "%s, no predicate, no path search: §7.6)" % (labels, mode)
    if phase == "selection" and not predicate:
        return "a selection stop without a predicate (§19.8)"
    if reason == "max_contexts":
        if phase == "selection" and (long_ or not predicate):
            return "a selection stop at max_contexts without a predicate on a pattern of L <= k (§19.6)"
        if phase == "discovery" and predicate and not long_ and \
                selection_pass in ("completed", "stopped", "not_admitted", "not_started"):
            return "a discovery stop at max_contexts with a predicate's pass (its raw threshold is " \
                   "max_predicate_contexts: §19.6)"
    if reason == "max_predicate_contexts":
        if phase == "discovery" and long_:
            return "a discovery stop at max_predicate_contexts of a pattern longer than k (§20.9: the held walks, " \
                   "phase extension)"
        if phase == "extension" and ls != "supported_paths":
            return "an extension stop at max_predicate_contexts without a supported-path search (§20.9)"
    if reason == "max_annotation_work" and phase in ("extension", "selection") and \
            not (long_ and ls == "supported_paths"):
        return "a %s stop at max_annotation_work without a supported-path search (§20.5, §20.9)" % phase
    return None


def _threshold_stop_count(reason: str, phase: Any, cnt: Callable[[str], Optional[dict]], main: Optional[dict],
                          predicate: bool, long_: bool) -> Optional[dict]:
    """The count a stop_at_threshold stop of |reason| interrupted: max_contexts the raw contexts in discovery (a
    constant normal form runs the engine's own threshold), the selected ones in the selection (§19.6);
    max_predicate_contexts the raw contexts in discovery, the supported paths whose walks it held (§20.9);
    max_anchors the anchors; max_paths the paths (with a predicate the selected supported paths, §20.9)."""
    if reason == "max_contexts":
        return cnt("selected") if phase == "selection" else cnt("contexts")
    if reason == "max_predicate_contexts":
        return cnt("contexts") if not long_ else (cnt("supported_paths") or cnt("paths"))
    if reason == "max_anchors":
        return cnt("anchors")
    if reason == "max_paths":
        return cnt("selected") if predicate and cnt("selected") else main
    return None


def _threshold_value_problem(sc: dict, body: dict, reason: str, phase: Any, stopped: Optional[dict],
                             predicate: bool, long_: bool) -> Optional[str]:
    """Why a discovery stop of stop_at_threshold on the contexts of a pattern of L <= k cannot have happened with the
    count it states, or None: on a masked BASIC or CANONICAL graph discovery stops once its running lower bound
    passes the threshold (max_contexts, with a predicate max_predicate_contexts: §7.5, §19.6), and the `at_least`
    value is that bound (the engine's threshold_crossed compares the sum of the counted ranges' lower bounds, which
    the halted count states: pattern_search.cpp running_lower, tally). Not judged on a graph without its mask (the
    running upper bound decides, the value is the lower parts only: fixture unmasked_stop_at_threshold, at_least 0)
    or on a PRIMARY graph (an even-k wrapped union's bound need not be the stated value: §18)."""
    idx = body.get("index") if isinstance(body.get("index"), dict) else {}
    if long_ or phase not in (None, "discovery") or idx.get("counting") == "upper_bound" or \
            idx.get("graph_mode") not in ("basic", "canonical"):
        return None
    if reason == "max_contexts" and not predicate:
        cap_name = "max_contexts"
    elif reason == "max_predicate_contexts" and predicate:
        cap_name = "max_predicate_contexts"
    else:
        return None
    cap = _number(_limit(sc, body, cap_name))
    if cap is None or not isinstance(stopped, dict) or stopped.get("relation") != "at_least" or \
            not isinstance(stopped.get("value"), (int, float)) or stopped["value"] > cap:
        return None
    return "a %s stop (stop_at_threshold) with contexts %s, not past %s %d (§7.5: the running lower bound passed " \
           "it)" % (reason, _count_text(stopped), cap_name, int(cap))


def _hiding(sc: dict, ans: dict, i: int, e: dict) -> Tuple[str, List[str]]:
    """(what hides part of entry i's list, what is wrong with the entry's stated reasons). Only a withheld list or a
    cut shortens a list of contexts or paths (SPEC-pattern-search §7.5: in partial every stop of the graph search
    states a cut, in all_or_count a list not handed over is withheld; a stop of the label reads or retrieval_complete
    false alone leave every context listed); '' when the list holds every result."""
    if e.get("retrieval_complete") is True:
        hidden = ""
    else:
        parts = []
        if e.get("withheld"):
            parts.append("withheld %s" % (e["withheld"].get("reason") if isinstance(e["withheld"], dict) else e["withheld"]))
        if e.get("cut"):
            parts.append("cut %s" % (e["cut"].get("reason") if isinstance(e["cut"], dict) else e["cut"]))
        hidden = ", ".join(parts)
    body = ans.get("body") if isinstance(ans.get("body"), dict) else {}
    return hidden, _pattern_reason_problems(sc, body, i, e)


def _resolve_reserve_problem(body: dict) -> Optional[str]:
    """Why a /resolve answer's stated finalisation reserve is not the host's, or None: `limits`, stated with every
    budget (SPEC-labeled-traversal-core §4.5: "the effective budget, the reserve, the clamps"), states
    finalize_reserve_ms, the host's 250 ms (fixed in this build; capabilities resolve.time_budget
    finalize_reserve_configurable false). An answer without `limits` (no budget) states none. Read by invariants on
    every answer and by the stop rule (_resolve_stop_problem), never believed in place of the host's."""
    lim = body.get("limits") if isinstance(body, dict) else None
    if not isinstance(lim, dict):
        return None
    if "finalize_reserve_ms" not in lim:
        return "limits states no finalize_reserve_ms (SPEC §4.5: the reserve is stated with the budget, the host's " \
               "%g ms)" % _RESOLVE_RESERVE_MS
    v = lim["finalize_reserve_ms"]
    if not _is_num(v) or v != _RESOLVE_RESERVE_MS:
        return "limits.finalize_reserve_ms %s, not the host's %g ms (SPEC §4.5: fixed in this build, not " \
               "configurable)" % (json.dumps(v), _RESOLVE_RESERVE_MS)
    return None


def _resolve_stop_problem(sc: dict, body: dict, client: Optional[dict] = None) -> Optional[str]:
    """Why a /resolve answer's stated stop cannot have happened, or None (SPEC-labeled-traversal-core §4.5: the
    work stops at the budget less the reserve, the only stop reason being time). The deadline starts once the body
    is parsed, timing.elapsed_ms's clock later (traverse.cpp process_resolve_request: after the request's fields
    and the LabelOracle), so elapsed_ms is no measure of the deadline's time: a valid stop can state far less than
    the budget less the reserve when that gap is long (a stall). Only a clock whose interval holds the deadline's
    can show a stop before the work deadline: the client's wall (|client| wall_ms: a direct runner's send to last
    byte, the MCP arm's service-measured backend exchange, mcp_phases.backend_ms). Without one the time is not
    judged. /traverse and /pattern read their stop and their elapsed time on one clock, so their rules keep it."""
    stop = body.get("stop") if isinstance(body, dict) else None
    if not isinstance(stop, dict):
        return None
    if stop.get("reason") != "time":
        return "a resolve stop of reason %r (SPEC §4.5 states time only)" % stop.get("reason")
    lim = body.get("limits") if isinstance(body.get("limits"), dict) else {}
    budget = _effective_budget(sc, body)
    if budget is None:
        return "a resolve time stop without a time budget"
    stated = _resolve_reserve_problem(body)
    if stated and "finalize_reserve_ms" in lim:
        return "a resolve time stop with %s" % stated
    reserve = _RESOLVE_RESERVE_MS
    wall = _number((client or {}).get("wall_ms"))
    if wall is not None and wall < budget - reserve - 1.0:
        return "a resolve time stop within a client wall of %.0f ms, before the budget %.0f ms less the reserve " \
               "%.0f ms (the client's interval holds the deadline's)" % (wall, budget, reserve)
    return None


# the cap trigger reasons of a walk and the knob each one's walk_domain limitation names (src/cli/traverse.cpp
# cap_knob); resource_limit names the knob of its resource_stop
_CAP_KNOBS = {"max_steps": "bounds.max_steps", "max_live_paths": "bounds.max_live_paths",
              "beam_pruned": "bounds.max_live_paths", "max_paths": "bounds.max_paths",
              "max_output_bp": "bounds.max_output_bp", "time_budget": "bounds.time_budget_ms"}
_RESOURCE_KNOBS = {"work": "bounds.max_work_units", "memory": "bounds.max_memory_mb", "cancelled": "attempt_id",
                   "attempt_deadline": "attempt_id"}
# a cap trips when the next head's followed successors (at most 4) would pass it (walker.cpp cap_check)
_MAX_FOLLOWED = 4


def _arm_stop_problem(sc: dict, ans: dict, r: dict, i: int, side: str) -> Optional[str]:
    """Why arm |side| of seed result |r| cannot have been stopped by its stated cap trigger, or None. The trigger's
    walk_domain limitation must name its knob with the observed demand past the limit (src/cli/traverse.cpp
    arm_limitations: the trigger's observed is what the cap compared), the arm's counters must stand at the cap
    (max_output_bp: the arm's output_bp; max_steps: the seed's steps; within the at most 4 successors a head
    follows), a time stop at or after the budget (walker.cpp time_exceeded), resource_limit with the seed's
    resource_stop. Judged on every answer shape: a graphlet result (`detail: graphlet`, and the MCP's export) carries
    its arms' cap_trigger, limitations and counters as JSON beside the MGT text (traverse.cpp arm_certificate_json:
    every detail level), so only a quantity the answer lacks is not judged (a missing limitations list, counter or
    timing)."""
    armo = (r.get("arms") or {}).get(side) or {}
    trig = (armo.get("cap_trigger") or {}).get("reason")
    if trig not in _RESOURCE_ENDS:
        return None
    body = _body(ans)
    req = sc.get("request") if isinstance(sc.get("request"), dict) else {}
    bounds = dict((req.get("strategy") or {}).get("bounds") or {})
    bounds.update(((body.get("strategy") or {}).get("bounds") or {}) if isinstance(body, dict) else {})
    rs = r.get("resource_stop") if isinstance(r.get("resource_stop"), dict) else None
    if trig == "resource_limit":
        if rs is None:
            return "%s: truncated by resource_limit without a resource_stop" % side
        knob = _RESOURCE_KNOBS.get(rs.get("resource"))
    else:
        knob = _CAP_KNOBS.get(trig)
    if trig == "time_budget":
        budget = _effective_budget(sc, body)
        if budget is not None and isinstance(r.get("timing"), dict):
            latest, _ = _traverse_time_stop(r, budget)
            if latest is not None and latest < budget - 1.0:
                return "%s: truncated by time_budget by %.0f ms, before the %.0f ms budget" % (side, latest, budget)
    lims = armo.get("limitations")
    if isinstance(lims, list) and knob:
        mine = [x for x in lims if x.get("kind") == "walk_domain" and x.get("knob") == knob]
        if not mine:
            return "%s: truncated by %s, but no walk_domain limitation names %s" % (side, trig, knob)
        obs, lim = _number(mine[0].get("observed")), _number(mine[0].get("limit"))
        if obs is not None and lim is not None and knob != "attempt_id":
            if (obs < lim - 1.0) if trig == "time_budget" else (obs <= lim):
                return "%s: truncated by %s with observed %s within the limit %s of %s" % (
                    side, trig, mine[0].get("observed"), mine[0].get("limit"), knob)
    if trig == "max_output_bp":
        out_bp, cap = _number((armo.get("counters") or {}).get("output_bp")), _number(bounds.get("max_output_bp"))
        if out_bp is not None and cap is not None and out_bp + _MAX_FOLLOWED <= cap:
            return "%s: truncated by max_output_bp at output_bp %.0f of max_output_bp %.0f" % (side, out_bp, cap)
    if trig == "max_steps":
        steps = [_number(((x or {}).get("counters") or {}).get("steps")) for x in (r.get("arms") or {}).values()]
        cap = _number(bounds.get("max_steps"))
        if steps and None not in steps and cap is not None and sum(steps) + _MAX_FOLLOWED <= cap:
            return "%s: truncated by max_steps after %.0f steps of max_steps %.0f" % (side, sum(steps), cap)
    if trig == "max_paths":
        cap = _number(bounds.get("max_paths"))
        n = len(armo.get("paths") or [])
        if cap is not None and isinstance(armo.get("paths"), list) and n + 2 * _MAX_FOLLOWED <= cap and \
                isinstance(armo.get("frontier_remaining"), dict) and \
                n + (_number(armo["frontier_remaining"].get("live_paths")) or 0) + _MAX_FOLLOWED <= cap:
            return "%s: truncated by max_paths with %d paths and %s live of max_paths %.0f" % (
                side, n, armo["frontier_remaining"].get("live_paths"), cap)
    return None


def _resource_end_problem(r: dict, armo: dict, side: str, reasons: Iterable[str]) -> Optional[str]:
    """Why a path of arm |side| that ended by the resource reason(s) |reasons| cannot have been ended by them, or
    None. A resource reason ends a path only on an arm a cap left short of its radius (SPEC-labeled-traversal-core
    §6.7, §7.0: `truncated`, or `pruned` by a beam, with `cap_trigger` naming the cap), and the path's reason is the
    cap in cap_trigger or a later cap that ended walks after it, which the arm's walk_domain limitations name (§7.0
    table: "the cap in cap_trigger, then any other cap that ended walks after it")."""
    reasons = set(reasons) & _RESOURCE_ENDS
    if not reasons:
        return None
    named = "/".join(sorted(reasons))
    status = armo.get("status")
    trig = (armo.get("cap_trigger") or {}).get("reason")
    if status not in ("truncated", "pruned"):
        return "%s: a path ended by %s on an arm %s" % (side, named, status or "without a status")
    if trig not in _RESOURCE_ENDS:
        return "%s: a path ended by %s on a %s arm without a resource cap trigger (cap_trigger %r)" % (
            side, named, status, trig)
    if trig in reasons:
        return None
    lims = armo.get("limitations")
    if not isinstance(lims, list):
        return None   # the arm states no limitations: a later cap cannot be told from a wrong reason
    knobs = {x.get("knob") for x in lims if isinstance(x, dict) and x.get("kind") == "walk_domain"}
    rs = r.get("resource_stop") if isinstance(r.get("resource_stop"), dict) else {}
    for reason in reasons:
        knob = _RESOURCE_KNOBS.get(rs.get("resource")) if reason == "resource_limit" else _CAP_KNOBS.get(reason)
        if knob and knob in knobs:
            return None
    return "%s: a path ended by %s, but the arm's cap trigger is %s and no walk_domain limitation names a later cap" % (
        side, named, trig)


def _resource_end_below_frontier(armo: dict, side: str, p: dict) -> Optional[str]:
    """Why path |p| cannot have ended by the resource reason it states, or None: a cap ends walks at the level it
    trips, at or past the arm's complete_to_bp, and "every path that ended before complete_to_bp ended for the
    reported semantic reason" (SPEC-labeled-traversal-core §6.10); not judged where the arm states no
    complete_to_bp."""
    ends = ({p.get("path_reason")} | set(p.get("end_reasons") or {})) & _RESOURCE_ENDS
    ctb, L = armo.get("complete_to_bp"), p.get("length_bp")
    if not ends or not _finite(ctb) or not _finite(L) or L >= ctb:
        return None
    return "%s: path %s ended by %s at %d bp, below the arm's complete_to_bp %d (every path that ended before it " \
           "ended for a semantic reason: §6.10)" % (side, p.get("id"), "/".join(sorted(ends)), L, ctb)


def _path_resource_problems(sc: dict, r: dict, i: int) -> List[str]:
    """_resource_end_problem for every path of the requested arms of seed result |r| (its path_reason and the label
    end reasons it counts: a resource reason censors the path's labels with the same reason, §7.0), and
    _resource_end_below_frontier."""
    probs: List[str] = []
    for side in _arms_of(sc, i):
        armo = (r.get("arms") or {}).get(side)
        if not isinstance(armo, dict):
            continue
        seen = set()
        for p in armo.get("paths") or []:
            ends = ({p.get("path_reason")} | set(p.get("end_reasons") or {})) & _RESOURCE_ENDS
            if ends and frozenset(ends) not in seen:
                seen.add(frozenset(ends))
                m = _resource_end_problem(r, armo, side, ends)
                if m:
                    probs.append("seed %d %s (path %s)" % (i, m, p.get("id")))
            m = _resource_end_below_frontier(armo, side, p)
            if m:
                probs.append("seed %d %s" % (i, m))
    return probs


def _traverse_reason_problems(sc: dict, ans: dict) -> List[str]:
    """Every requested arm's stated cap trigger that cannot have happened (_arm_stop_problem), over the walked seeds."""
    probs: List[str] = []
    try:
        res = _get(_body(ans), "results")
    except Missing:
        return probs
    for i in range(len(res)):
        r = _result(ans, i)
        if "arms" not in r:
            continue
        for side in _arms_of(sc, i):
            m = _arm_stop_problem(sc, ans, r, i, side)
            if m:
                probs.append("seed %d %s" % (i, m))
    return probs


def time_limited(ans: dict, sc: Optional[dict] = None, client: Optional[dict] = None, arm: Optional[str] = None) -> str:
    """What makes an answer time-limited (SPEC-labeled-traversal-core §6.8: determinism holds provided no time
    budget tripped; SPEC-pattern-search §7.9: determinism "time_limited"); '' when nothing did. With the scenario
    |sc| (and the runner's |client| and |arm|), also the front ingress's 504 at its cut (ingress_cut_allowed): the
    budget ran out before the host's answer could be seen, a time limitation as the 503 deadline is (VERIFY-V6-2
    L1). The runners record it beside an execution's digest and pass {digest, time_limited} for the earlier repeats
    (`previous`)."""
    own = _time_limited(ans)
    if not own and sc is not None and ingress_cut_allowed(sc, ans, client, arm):
        return INGRESS_CUT_LIMITATION
    return own


def _time_limitation_problem(sc: dict, ans: dict, arm: Optional[str] = None,
                             client: Optional[dict] = None) -> Optional[str]:
    """Why an answer's stated time limitation cannot have happened, or None: a /pattern time stop or limitation
    before its work time, a /resolve stop before its budget less the reserve by the client's wall (|client|,
    _resolve_stop_problem), a /traverse time stop before its
    budget, an attempt-bound stop before its bound (_attempt_stops; AttemptAtBound in a 200 seed result only on the
    MCP arm, which relays the host's 503 so)."""
    body = ans.get("body")
    if not isinstance(body, dict):
        return None
    route = sc.get("route", "")
    if route == "POST /pattern":
        for i, e in enumerate(body.get("patterns") or []):
            if isinstance(e, dict) and (e.get("determinism") == "time_limited" or (e.get("stop") or {}).get("reason") == "time"):
                return _pattern_time_problem(sc, body, "p%d: a time limitation" % i)
        return None
    if route == "POST /resolve":
        return _resolve_stop_problem(sc, body, client)
    if route == "POST /traverse":
        probs = [m for m in _traverse_reason_problems(sc, ans) if "time_budget" in m]
        probs += [p for _, p in _attempt_stops(sc, ans, "mcp").values() if p]
        if arm not in (None, "mcp"):
            probs += ["seed %d states the host's AttemptAtBound in a 200 (a direct arm receives it as a 503)" % i
                      for i, r in enumerate(body.get("results") or [])
                      if isinstance(r, dict) and "arms" not in r and _AT_BOUND_RE.fullmatch(_seed_error_text(r) or "")]
        return probs[0] if probs else None
    return None


# ---------------------------------------------------------------------------------------------- pattern checks
def check_pattern_slot_error(sc, ans, a, arm):
    e = _entry(ans, a["pattern"])
    err = e.get("error")
    if not isinstance(err, dict):
        return FAIL, "entry answered (no error slot), expected %s" % a["code"]
    if err.get("code") != a["code"]:
        return FAIL, "slot code %r, expected %r (%s)" % (err.get("code"), a["code"], err.get("message"))
    msg = err.get("message") or ""
    if "text" in a and msg != a["text"]:
        return FAIL, "slot message %r, expected %r" % (msg, a["text"])
    if "text_prefix" in a and not msg.startswith(a["text_prefix"]):
        return FAIL, "slot message %r does not start with %r" % (msg, a["text_prefix"])
    for frag in a.get("text_contains", []):
        if frag not in msg:
            return FAIL, "slot message %r lacks %r" % (msg, frag)
    for path, want in (a.get("fields") or {}).items():
        got = _get(e, path)
        if got != want:
            return FAIL, "slot %s = %r, expected %r" % (path, got, want)
    if "counts" in e or "results" in e:
        return FAIL, "a refused slot carries counts or results"
    return PASS, "slot %s" % a["code"]


# the graph counts (SPEC-pattern-search §7.4, §7.7, §12.1, §20.4): at_least only after a stop of the graph search; a
# stop of the labels' reads or of the output leaves them as the search made them (§7.5: "the graph count is then still
# the discovery's")
_GRAPH_COUNTS = ("contexts", "anchors", "paths", "supported_paths")
_LABEL_PHASES = {"label_discovery", "placement", "output"}


def _stopless_inexact(body: dict, e: dict, req: dict, name: str, rel: str) -> Optional[str]:
    """Why graph count |name| of entry |e| may have the inexact relation |rel| with no stop that leaves it so, by the
    SPEC, or None: bounds on a graph without its mask (SPEC-pattern-search §7.4); on a supported-path search,
    supported_paths and paths at_least after a predicate's pruning and paths at_least after any branch pruned before
    the pattern's end (§20.4, §20.9). Otherwise a graph count is at_least or bounds only after a stop (§7.4, §12.1)."""
    counts = e.get("counts") or {}
    if rel == "bounds" and ((body.get("index") or {}).get("counting") == "upper_bound"):
        return "a graph without its mask (index.counting upper_bound, SPEC-pattern-search §7.4)"
    if name in ("supported_paths", "paths") and req.get("long_search") == "supported_paths" and rel == "at_least":
        sp = counts.get("supported_paths") if isinstance(counts.get("supported_paths"), dict) else {}
        by_predicate = _number(sp.get("branches_pruned_by_predicate")) or 0
        pruned = _number(sp.get("branches_pruned")) or 0
        if by_predicate > 0:
            return "the predicate pruned %d branch(es) of the search (SPEC-pattern-search §20.4, §20.9)" % by_predicate
        if name == "paths" and pruned > 0:
            return "%d branch(es) pruned before the pattern's end (SPEC-pattern-search §20.4)" % pruned
    return None


def check_count_relation(sc, ans, a, arm):
    """A count's relation and value (`relation`, `value`, `min`, `max`; `level` for by_mismatches). A count a stop
    left inexact (at_least, unknown or bounds beside a stated stop of a reason that can leave it so,
    _INEXACT_STOP_REASONS; for a graph count not a stop of the labels' reads or the output, which leave it as the
    search made it), or a graph count the SPEC states inexact without such a stop (_stopless_inexact), cannot show an
    exact expectation: inconclusive, unless the scenario expects that relation or the bound it states decides
    (at_least v >= min passes, at_least v > max fails, a value outside the bounds fails). A graph count at_least
    below a minimum with nothing that can leave it so contradicts the SPEC (§7.4: at_least means a stop interrupted
    it) and fails; another count (labels, occurrences, tested, selected) has relation rules of its own (§12.1,
    §14.6, §19.7): at_least below a minimum is inconclusive there."""
    e = _entry(ans, a.get("pattern", 0))
    if "error" in e:
        return FAIL, "entry refused (%s)" % e["error"].get("code")
    name = a.get("count", "contexts")
    c = _get(e, "counts." + name)
    if "level" in a:
        c = _get(c, "by_mismatches.%d" % a["level"])
    value, rel = _count(c)
    stop = e.get("stop") or {}
    graph = name in _GRAPH_COUNTS
    inexact = rel in ("at_least", "unknown", "bounds")
    stopped = inexact and stop.get("reason") in _INEXACT_STOP_REASONS and not (
        graph and stop.get("phase") in _LABEL_PHASES)
    if stopped:
        bad = _pattern_reason_problems(sc, _body(ans), a.get("pattern", 0), e)
        if bad:
            return FAIL, "%s %s %s after a %s/%s stop that cannot have happened: %s" % (
                name, rel, value, stop.get("phase"), stop.get("reason"), bad[0])
    req = sc.get("request") if isinstance(sc.get("request"), dict) else {}
    body = _body(ans)
    cause = ("the stop %s/%s" % (stop.get("phase"), stop.get("reason"))) if stopped else \
        (_stopless_inexact(body if isinstance(body, dict) else {}, e, req, name, rel) if inexact and graph else None)
    if "relation" in a:
        rels = a["relation"] if isinstance(a["relation"], list) else [a["relation"]]
        if rel not in rels:
            if cause:
                return INCONCLUSIVE, "%s %s %s after %s: not the relation expected (%s)" % (
                    name, rel, value, cause, rels)
            return FAIL, "%s relation %s (value %s), expected %s" % (name, rel, value, rels)
    lower = c.get("lower", value) if rel == "bounds" else value
    upper = c.get("upper") if rel == "bounds" else value if rel == "exact" else None
    if "value" in a:
        if rel == "exact" and value != a["value"]:
            return FAIL, "%s = exact %s, expected exact %s" % (name, value, a["value"])
        if rel != "exact":
            if (lower is not None and lower > a["value"]) or (upper is not None and upper < a["value"]):
                return FAIL, "%s = %s %s%s, which excludes the exact %s expected" % (
                    name, rel, value, " [%s, %s]" % (lower, upper) if rel == "bounds" else "", a["value"])
            if cause:
                return INCONCLUSIVE, "%s %s %s after %s: the exact %s expected not shown" % (
                    name, rel, value, cause, a["value"])
            return FAIL, "%s = %s %s, expected exact %s" % (name, rel, value, a["value"])
    if "min" in a:
        if rel == "unknown":
            if cause:
                return INCONCLUSIVE, "%s unknown after %s" % (name, cause)
            return FAIL, "%s unknown, expected at least %d" % (name, a["min"])
        if rel == "at_least" and value is not None and value < a["min"]:
            if cause:
                return INCONCLUSIVE, "%s at least %d after %s, %d expected" % (name, value, cause, a["min"])
            if not graph:
                return INCONCLUSIVE, "%s at least %d, %d expected: a count of relation rules of its own " \
                                     "(SPEC-pattern-search §12.1, §14.6, §19.7)" % (name, value, a["min"])
            return FAIL, "%s at least %d below the minimum %d, with no stop that leaves the count inexact " \
                         "(SPEC-pattern-search §7.4: at_least means a stop interrupted it)" % (name, value, a["min"])
        if rel == "exact" and value < a["min"]:
            return FAIL, "%s exact %d below %d" % (name, value, a["min"])
        if rel == "bounds":
            if c.get("upper", 0) < a["min"]:
                return FAIL, "%s bounds upper %s below %d" % (name, c.get("upper"), a["min"])
            if lower is not None and lower < a["min"]:
                return INCONCLUSIVE, "%s bounds [%s, %s] straddle the minimum %d%s" % (
                    name, lower, c.get("upper"), a["min"], " (after %s)" % cause if cause else "")
    if "max" in a and rel in ("exact", "at_least", "bounds") and lower is not None and lower > a["max"]:
        return FAIL, "%s %s %s above %d" % (name, rel, lower, a["max"])
    return PASS, "%s %s %s" % (name, rel, value)


def check_withheld(sc, ans, a, arm):
    e = _entry(ans, a.get("pattern", 0))
    w = e.get("withheld", "<absent>")
    if w == "<absent>":
        raise Missing("withheld")
    reason = w.get("reason") if isinstance(w, dict) else None
    want = a.get("reason")
    wants = want if isinstance(want, list) else [want]
    if reason not in wants:
        return FAIL, "withheld %r, expected %r" % (reason, want)
    return PASS, "withheld %r" % reason


def check_cut(sc, ans, a, arm):
    e = _entry(ans, a.get("pattern", 0))
    c = e.get("cut", "<absent>")
    if c == "<absent>":
        raise Missing("cut")
    reason = c.get("reason") if isinstance(c, dict) else None
    want = a.get("reason")
    wants = want if isinstance(want, list) else [want]
    if reason not in wants:
        return FAIL, "cut %r, expected %r" % (reason, want)
    return PASS, "cut %r" % reason


def check_stop(sc, ans, a, arm):
    e = _entry(ans, a.get("pattern", 0))
    if "stop" not in e:
        raise Missing("stop")
    s = e["stop"]
    if a.get("null"):
        return (PASS, "stop null") if s is None else (FAIL, "stop %r, expected null" % s)
    if s is None:
        return FAIL, "stop null, expected %s/%s" % (a.get("phase"), a.get("reason"))
    for f in ("phase", "reason"):
        if f in a:
            wants = a[f] if isinstance(a[f], list) else [a[f]]
            if s.get(f) not in wants:
                return FAIL, "stop.%s %r, expected %r" % (f, s.get(f), a[f])
    if "mismatch_level_in" in a and s.get("mismatch_level") not in a["mismatch_level_in"]:
        return FAIL, "stop.mismatch_level %r, expected one of %r" % (s.get("mismatch_level"), a["mismatch_level_in"])
    return PASS, "stop %s/%s" % (s.get("phase"), s.get("reason"))


def check_retrieval_complete(sc, ans, a, arm):
    e = _entry(ans, a.get("pattern", 0))
    v = e.get("retrieval_complete", "<absent>")
    if v == "<absent>":
        raise Missing("retrieval_complete")
    if v != a["value"]:
        return FAIL, "retrieval_complete %r, expected %r" % (v, a["value"])
    return PASS, "retrieval_complete %r" % v


def check_notes(sc, ans, a, arm):
    e = _entry(ans, a.get("pattern", 0))
    notes = e.get("notes")
    if notes is None:
        raise Missing("notes")
    for n in a.get("has", []):
        if n not in notes:
            return FAIL, "note %r missing (notes %r)" % (n, notes)
    for n in a.get("lacks", []):
        if n in notes:
            return FAIL, "note %r present" % n
    return PASS, "notes %r" % notes


_SEARCH_BUDGET_STOPS = {"max_steps", "time", "max_annotation_work", "max_memory"}
# the engine's phases (SPEC-pattern-search §7.6: the range search, the deferred scans, the release and its delivery,
# a path search's extension), whose stop the engine's low-complexity diagnostic follows; the labels' reads, a
# predicate's selection and the output come after it (§7.8)
_ENGINE_PHASES = {"discovery", "mask_scan", "extraction", "extension"}


def _note_stop(sc: dict, e: dict) -> str:
    """How entry |e|'s stop bears on the low_complexity_pattern note (SPEC-pattern-search §7.8, §8.10, §20.5): "budget"
    for a budget stop of the search (max_steps, time, the supported-path search's max_annotation_work and max_memory)
    in an engine phase, a stop without a phase counted as the engine's, which leaves the note out; "selection" for a
    selection-phase stop of a supported-path search, whose selection runs in the search (a held walk the account
    refused) or after it, so that the note may be either; "" for none, a stop_at_threshold stop or a later phase's
    stop (labels, selection, output), which leave the note as it is."""
    stop = e.get("stop") if isinstance(e.get("stop"), dict) else {}
    req = sc.get("request") if isinstance(sc.get("request"), dict) else {}
    phase = stop.get("phase")
    if stop.get("reason") in _SEARCH_BUDGET_STOPS and (phase is None or phase in _ENGINE_PHASES):
        return "budget"
    if phase == "selection" and req.get("long_search") == "supported_paths":
        return "selection"
    return ""


def check_low_complexity_note(sc, ans, a, arm):
    """The `low_complexity_pattern` note as SPEC-pattern-search §7.8 states it: on an exact pattern that sdust flags
    (args.expected, computed by the builder with seqtools.sdust, the host's T = 20, W = 64), when the answer states
    no budget stop of the search (_note_stop: an engine phase's); never beside such a stop. A stop of a later phase
    (the labels' reads, a predicate's selection, the output) leaves the note as the engine decided it, so the note
    follows sdust there; a supported-path search's selection stop allows either, a present note only where sdust
    flags the pattern. A pattern of more than 191 bases whose diagnostic the work time cut (determinism
    time_limited) may lack it: inconclusive then, as is an expected note a budget stop left out."""
    e = _entry(ans, a.get("pattern", 0))
    if "error" in e:
        return FAIL, "entry refused (%s)" % e["error"].get("code")
    notes = e.get("notes")
    if notes is None:
        raise Missing("notes")
    has = "low_complexity_pattern" in notes
    stop = e.get("stop") or {}
    kind = _note_stop(sc, e)
    if not has and (stop.get("reason") in _SEARCH_BUDGET_STOPS or e.get("determinism") == "time_limited"):
        bad = _pattern_reason_problems(sc, _body(ans), a.get("pattern", 0), e)
        if bad:
            return FAIL, "no note after a stop or time limitation that cannot have happened: %s" % bad[0]
    if kind == "budget":
        if has:
            return FAIL, "low_complexity_pattern beside a budget stop %s/%s" % (stop.get("phase"), stop.get("reason"))
        if a["expected"]:
            return INCONCLUSIVE, "a budget stop %s/%s leaves the note out" % (stop.get("phase"), stop.get("reason"))
        return PASS, "no note (a budget stop, and sdust does not flag the pattern)"
    if kind == "selection" and not has and a["expected"]:
        return INCONCLUSIVE, "a supported-path search's selection stop %s/%s, which can come from the search's own " \
                             "refusal of a held walk: the note may be left out" % (
                                 stop.get("phase"), stop.get("reason"))
    later = " beside a later phase's stop %s/%s" % (stop.get("phase"), stop.get("reason")) if \
        stop.get("phase") is not None and stop.get("phase") not in _ENGINE_PHASES else ""
    if has == bool(a["expected"]):
        return PASS, "note %s, as sdust says%s" % ("present" if has else "absent", later)
    try:
        L = _pattern_length(*_request_pattern(sc, a.get("pattern", 0)))
    except Missing:
        L = 0
    if not has and L > 191 and e.get("determinism") == "time_limited":
        return INCONCLUSIVE, "the diagnostic of a %d-base pattern was cut by the work time" % L
    return FAIL, "low_complexity_pattern %s, sdust %s the pattern%s" % (
        "present" if has else "absent", "flags" if a["expected"] else "does not flag", later)


def check_stops_stated(sc, ans, a, arm):
    """Every answered entry either kept the budget (stop null, its primary count exact) or states a stop of the
    reasons args.reasons (a time budget: "time"); a discovery or mask-scan stop leaves the count inexact
    (SPEC-pattern-search §7.6: at_least or unknown, the stop sticky for the later patterns, which
    counts_consistent checks)."""
    pats = _get(_body(ans), "patterns")
    n_stop = 0
    for i, e in enumerate(pats):
        if "error" in e or "counts" not in e:
            continue
        stop = e.get("stop")
        name = _primary_count_name(e, sc["request"])
        c = e["counts"].get(name) or {}
        if stop is None:
            if c.get("relation") != "exact":
                return FAIL, "p%d: no stop, but %s %s" % (i, name, c.get("relation"))
            continue
        n_stop += 1
        if stop.get("reason") not in a["reasons"]:
            return FAIL, "p%d: stop %s/%s, expected one of %s" % (i, stop.get("phase"), stop.get("reason"), a["reasons"])
        if stop.get("phase") in ("discovery", "mask_scan") and c.get("relation") == "exact":
            return FAIL, "p%d: an exact %s after a %s stop" % (i, name, stop.get("phase"))
    return PASS, "%d of %d entries stopped by %s, the others exact" % (n_stop, len(pats), a["reasons"])


def check_results_levels(sc, ans, a, arm):
    """Every released result lies at a level of at most args.max_level (mismatches first: a list cut inside the
    lowest level holds that level only), and at least args.min_returned results are released."""
    e = _entry(ans, a.get("pattern", 0))
    res = e.get("results")
    if res is None:
        raise Missing("results")
    if len(res) < a.get("min_returned", 1):
        return FAIL, "%d results released, expected at least %d" % (len(res), a.get("min_returned", 1))
    over = [r.get("mismatches", 0) for r in res if r.get("mismatches", 0) > a["max_level"]]
    if over:
        return FAIL, "%d of %d results above level %d (levels %s)" % (len(over), len(res), a["max_level"], sorted(set(over)))
    return PASS, "%d results, all at level <= %d" % (len(res), a["max_level"])


def check_absent(sc, ans, a, arm):
    """Exact 0 of the scope's count (contexts, or for a pattern longer than k its paths / anchors), and in a
    retrieval mode a complete empty list (SPEC §9: an exact 0 is an absence claim)."""
    i = a.get("pattern", 0)
    e = _entry(ans, i)
    if "error" in e:
        return FAIL, "entry refused (%s)" % e["error"].get("code")
    name = a.get("count") or ("contexts" if "contexts" in e.get("counts", {}) else "paths")
    value, rel = _count(_get(e, "counts." + name))
    if rel != "exact" or value != 0:
        return FAIL, "%s %s %s, expected exact 0" % (name, rel, value)
    if e.get("mode", sc["request"].get("mode")) != "count":
        if e.get("retrieval_complete") is not True:
            return FAIL, "exact 0 but retrieval_complete %r" % e.get("retrieval_complete")
        if e.get("results") not in ([], None) and len(e["results"]) != 0:
            return FAIL, "exact 0 with %d results" % len(e["results"])
    return PASS, "%s exact 0" % name


def _find_source(e: dict, t: dict, strand: str, level: Optional[int] = None,
                 scope: str = "any_offset") -> Tuple[List[Tuple[str, int]], List[dict]]:
    """The record's contexts of the truth window (in the request's scope: suffix keeps offset k - L only) and
    those of them the results hold, at |level| when given."""
    wanted = source_contexts(t)
    if scope == "suffix":
        wanted = [(km, p) for km, p in wanted if p == K - t["length"]]
    results = e.get("results")
    if results is None:
        raise Missing("results")
    found = []
    for kmer, p in wanted:
        for r in results:
            if r.get("kmer") == kmer and r.get("offset") == p and r.get("strand", r.get("orientation")) in (
                    [strand] if strand != "=" else ["=", "palindromic"]) and (
                    level is None or r.get("mismatches", 0) == level):
                found.append(r)
                break
    return wanted, found


def _hidden_outcome(sc: dict, ans: dict, i: int, e: dict, hidden_fmt: str, complete: str) -> Tuple[str, str]:
    """The outcome of a check whose result is missing from entry i's list: INCONCLUSIVE when a legitimate withheld
    list or cut hides part of it (hidden_fmt % what), FAIL when the stated reason cannot have happened or when the
    list holds every result (complete)."""
    hidden, bad = _hiding(sc, ans, i, e)
    if bad:
        return FAIL, (hidden_fmt % (hidden or "nothing hidden")) + "; a stated reason that cannot have happened: " + bad[0]
    if hidden:
        return INCONCLUSIVE, hidden_fmt % hidden
    return FAIL, complete


def check_contains_source_context(sc, ans, a, arm):
    i = a.get("pattern", 0)
    e = _entry(ans, i)
    if "error" in e:
        return FAIL, "entry refused (%s)" % e["error"].get("code")
    t = sc["truth"]
    if "truths" in sc and "truth" in a:
        t = sc["truths"][a["truth"]]
    strand = a.get("strand", t.get("strand", "+"))
    level = a.get("level")
    wanted, found = _find_source(e, t, strand, level, sc["request"].get("scope", "any_offset"))
    if not wanted:
        return FAIL, "the truth has no record k-mer holding the window (builder error)"
    require = a.get("require", "all")
    if require == "all" and len(found) != len(wanted):
        missing = [p for (km, p) in wanted if not any(r["kmer"] == km and r["offset"] == p for r in found)]
        detail = "%d of %d source contexts (strand %s, level %s) returned; missing offsets %s (returned %d" % (
            len(found), len(wanted), strand, level, missing[:8], len(e.get("results", [])))
        return _hidden_outcome(sc, ans, i, e, detail + "; %s)", detail + ", complete)")
    if require == "any" and not found:
        detail = "none of %d source contexts returned (strand %s, level %s" % (len(wanted), strand, level)
        return _hidden_outcome(sc, ans, i, e, detail + "; %s)", detail + ")")
    if "substitutions" in a:
        for r in found:
            if r.get("substitutions") != a["substitutions"]:
                return FAIL, "source context substitutions %r, expected %r" % (r.get("substitutions"), a["substitutions"])
    return PASS, "%d of %d source contexts returned (strand %s, level %s)" % (len(found), len(wanted), strand, level)


def check_source_absent(sc, ans, a, arm):
    """The source contexts are NOT among the results at any level (near misses): needs a complete list; nothing
    is assumed about the other contexts."""
    i = a.get("pattern", 0)
    e = _entry(ans, i)
    if "error" in e:
        return FAIL, "entry refused (%s)" % e["error"].get("code")
    t = sc["truth"]
    strand = a.get("strand", t.get("strand", "+"))
    wanted, found = _find_source(e, t, strand, None, sc["request"].get("scope", "any_offset"))
    if found:
        return FAIL, "%d source context(s) returned (level %s), expected out of reach" % (
            len(found), [r.get("mismatches") for r in found])
    hidden, bad = _hiding(sc, ans, i, e)
    if bad:
        return FAIL, "the source not among %d results; a stated reason that cannot have happened: %s" % (
            len(e.get("results") or []), bad[0])
    if hidden:
        return INCONCLUSIVE, "the source not among %d results, but the list is not complete (%s)" % (
            len(e.get("results") or []), hidden)
    return PASS, "none of %d source contexts among %d results (every context listed)" % (
        len(wanted), len(e.get("results", [])))


def _label_of(r: dict, column: str) -> Optional[dict]:
    for lab in r.get("labels") or []:
        if lab.get("column") == column:
            return lab
    return None


def check_source_label(sc, ans, a, arm):
    i = a.get("pattern", 0)
    e = _entry(ans, i)
    t = sc["truth"]
    strand = a.get("strand", t.get("strand", "+"))
    wanted, found = _find_source(e, t, strand, None, sc["request"].get("scope", "any_offset"))
    if not found:
        return _hidden_outcome(sc, ans, i, e, "no source context returned (%s)", "no source context returned")
    for r in found:
        if r.get("labels") is None:
            raise Missing("results[].labels")
        if _label_of(r, t["taxid"]):
            return PASS, "taxid %s among the labels of the source context (offset %d)" % (t["taxid"], r["offset"])
    return FAIL, "taxid %s on none of %d source contexts (labels_status %s)" % (
        t["taxid"], len(found), [r.get("labels_status") for r in found])


def check_source_placement(sc, ans, a, arm):
    """The record is placed: the source context's label (the record's taxid) lists an occurrence of the record at
    the window's 1-based coordinates and strand (SPEC §14.3). An occurrence_list the label does not state is a
    field the answer lacks (not exposed on the MCP arm); a null one states "not placed" (SPEC-pattern-search §14.5,
    the answer table), which fails on every arm where no other source context places the record (SD-377)."""
    i = a.get("pattern", 0)
    e = _entry(ans, i)
    t = sc["truth"]
    strand = a.get("strand", t.get("strand", "+"))
    wanted, found = _find_source(e, t, strand, None, sc["request"].get("scope", "any_offset"))
    if not found:
        return _hidden_outcome(sc, ans, i, e, "no source context returned (%s)", "no source context returned")
    coords = "%d-%d" % (t["start"] + 1, t["start"] + t["length"])
    unplaced = []   # the offsets of the source contexts whose label states occurrence_list null
    for r in found:
        lab = _label_of(r, t["taxid"])
        if lab is None:
            continue
        if "occurrence_list" not in lab:
            raise Missing("occurrence_list")
        occ = lab["occurrence_list"]
        if occ is None:
            unplaced.append(r.get("offset"))
            continue
        for o in occ:
            if o.get("record") == t["record"] and o.get("nt_coords") == coords and o.get("strand") == strand:
                m = _occurrence_problem(o, t)
                if m:
                    return FAIL, m
                return PASS, "%s %s %s placed (offset %d)" % (t["record"], coords, strand, r["offset"])
    if unplaced:
        return FAIL, "taxid %s states occurrence_list null (not placed, SPEC-pattern-search §14.5) on %d of %d " \
                     "source contexts (offsets %s): no occurrence %s %s %s" % (
                         t["taxid"], len(unplaced), len(found), unplaced[:6], t["record"], coords, strand)
    return FAIL, "no occurrence %s %s %s under taxid %s on %d source contexts" % (
        t["record"], coords, strand, t["taxid"], len(found))


def _coords_problem(o: dict) -> Optional[str]:
    """Why occurrence |o|'s 1-based inclusive nt_coords do not lie in its own record (1 <= start <= end <=
    nt_length, SPEC-pattern-search §14.3), or None (also where it states no coordinates or no length)."""
    m = re.fullmatch(r"(\d+)-(\d+)", str(o.get("nt_coords")))
    n = o.get("nt_length")
    if not m or not isinstance(n, int) or isinstance(n, bool):
        return None
    a, b = int(m.group(1)), int(m.group(2))
    if 1 <= a <= b <= n:
        return None
    return "occurrence of %s at %s on %s, outside its record of nt_length %d" % (o.get("record"), o.get("nt_coords"),
                                                                                  o.get("strand"), n)


def _occurrence_problem(o: dict, t: dict) -> Optional[str]:
    """What is wrong with the truth record's occurrence |o| (its record, coordinates and strand matched): its seq_id
    must be the index's for the record and its nt_length the record's length, which holds its coordinates."""
    if "seq_id" in o and o["seq_id"] != t.get("seq_id"):
        return "record %s placed with seq_id %s, the index's is %s" % (t["record"], o["seq_id"], t.get("seq_id"))
    if o.get("nt_length") not in (None, t.get("record_length")):
        return "nt_length %s, the record's length is %s" % (o.get("nt_length"), t.get("record_length"))
    return _coords_problem(o)


def _paths_with_sequence(e: dict, seq: str) -> List[dict]:
    res = e.get("results")
    if res is None:
        raise Missing("results")
    return [r for r in res if r.get("sequence") == seq]


def check_path_spells(sc, ans, a, arm):
    i = a.get("pattern", 0)
    e = _entry(ans, i)
    t = sc["truth"]
    strand = a.get("strand", t.get("strand", "+"))
    hits = [r for r in _paths_with_sequence(e, t["window"]) if r.get("strand", r.get("orientation")) in (
        [strand] if strand != "=" else ["=", "palindromic"])]
    if not hits:
        detail = "no path spells the record's bases (strand %s; %d paths returned, " % (strand, len(e.get("results", [])))
        return _hidden_outcome(sc, ans, i, e, detail + "%s)", detail + "every path listed)")
    for r in hits:
        if r.get("instance") not in (None, r["sequence"]) or r.get("offset", 0) != 0:
            return FAIL, "path instance/offset inconsistent"
        n = t["length"] - K + 1
        if "nodes" in r and len(r["nodes"]) != n:
            return FAIL, "path of %d nodes, expected %d" % (len(r["nodes"]), n)
    return PASS, "%d path(s) spell the record's bases" % len(hits)


def check_path_record_verified(sc, ans, a, arm):
    i = a.get("pattern", 0)
    e = _entry(ans, i)
    t = sc["truth"]
    strand = a.get("strand", t.get("strand", "+"))
    coords = "%d-%d" % (t["start"] + 1, t["start"] + t["length"])
    hits = _paths_with_sequence(e, t["window"])
    if not hits:
        detail = "no path spells the record's bases (%d paths, " % len(e.get("results", []))
        return _hidden_outcome(sc, ans, i, e, detail + "%s)", detail + "every path listed)")
    for r in hits:
        if r.get("labels") is None:
            continue
        lab = _label_of(r, t["taxid"])
        if lab is None:
            continue
        if lab.get("support") != "record_verified":
            return FAIL, "taxid %s carries the path with support %r" % (t["taxid"], lab.get("support"))
        # every occurrence the label verifies lies in its record; the record's own is placed with the index's
        # seq_id and the record's length (as source_placement places a context)
        outside = next((m for m in (_coords_problem(o) for o in lab.get("occurrence_list") or []) if m), None)
        if outside:
            return FAIL, "taxid %s: %s" % (t["taxid"], outside)
        for o in lab.get("occurrence_list") or []:
            if o.get("record") == t["record"] and o.get("nt_coords") == coords and o.get("strand") == strand:
                m = _occurrence_problem(o, t)
                if m:
                    return FAIL, m
                return PASS, "%s %s %s record_verified" % (t["record"], coords, strand)
        return FAIL, "taxid %s record_verified, but no occurrence %s %s %s (%s)" % (
            t["taxid"], t["record"], coords, strand, (lab.get("occurrence_list") or [])[:3])
    if all(r.get("labels") is None for r in hits):
        raise Missing("path labels")
    return FAIL, "taxid %s on none of the %d path(s) spelling the record's bases" % (t["taxid"], len(hits))


def check_path_absent(sc, ans, a, arm):
    i = a.get("pattern", 0)
    e = _entry(ans, i)
    seq = a.get("sequence") or sc["truth"].get("spliced_dna")
    hits = _paths_with_sequence(e, seq)
    if hits:
        return FAIL, "%d path(s) spell the sequence expected absent" % len(hits)
    hidden, bad = _hiding(sc, ans, i, e)
    if bad:
        return FAIL, "no path spells it among %d; a stated reason that cannot have happened: %s" % (
            len(e.get("results") or []), bad[0])
    if hidden:
        return INCONCLUSIVE, "no path spells it among %d, but the list is not complete (%s)" % (
            len(e.get("results") or []), hidden)
    return PASS, "no path spells it (every path listed: %d)" % len(e.get("results", []))


def check_predicate_selected(sc, ans, a, arm):
    """The source contexts are selected (in the results) or not (absent from a complete filtered list)."""
    i = a.get("pattern", 0)
    e = _entry(ans, i)
    if "error" in e:
        return FAIL, "entry refused (%s)" % e["error"].get("code")
    t = sc["truth"]
    strand = a.get("strand", t.get("strand", "+"))
    wanted, found = _find_source(e, t, strand, None, sc["request"].get("scope", "any_offset"))
    if a["selected"]:
        if len(found) != len(wanted):
            detail = "%d of %d source contexts selected" % (len(found), len(wanted))
            return _hidden_outcome(sc, ans, i, e, detail + " (%s)", detail)
        if "selection_strand_in" in a:
            for r in found:
                ss = r.get("selection_strands")
                if ss is None:
                    raise Missing("selection_strands")
                labs = r.get("selection_labels") or []
                col = a.get("label", t["taxid"])
                if col not in labs:
                    return FAIL, "label %s not in selection_labels %r" % (col, labs)
                if ss[labs.index(col)] not in a["selection_strand_in"]:
                    return FAIL, "selection strand %r for %s, expected %r" % (ss[labs.index(col)], col, a["selection_strand_in"])
        return PASS, "%d source contexts selected" % len(found)
    if found:
        return FAIL, "%d source context(s) selected, expected none" % len(found)
    hidden, bad = _hiding(sc, ans, i, e)
    if bad:
        return FAIL, "no source context selected; a stated reason that cannot have happened: %s" % bad[0]
    if hidden:
        return INCONCLUSIVE, "no source context selected, but the filtered list is not complete (%s)" % hidden
    return PASS, "source contexts not selected (every selected context listed)"


def check_predicate_either_mirror(sc, ans, a, arm):
    """Under predicate_strands "either": every returned '-' context whose reverse complement is a record k-mer
    of the truth is selected by the label with selection strand reverse_complement or both (SPEC §19.5), and at
    least args.min (default 1) such contexts are returned (the builder picks windows of inverted copies, whose
    mirrors the record holds); INCONCLUSIVE below that when a legitimate withheld list or cut hides them."""
    i = a.get("pattern", 0)
    e = _entry(ans, i)
    if "error" in e:
        return FAIL, "entry refused (%s)" % e["error"].get("code")
    t = sc["truth"]
    rec_kmers = {km for km, _ in source_contexts(t)}
    col = a.get("label", t["taxid"])
    checked = 0
    for r in e.get("results") or []:
        if r.get("strand") == "-" and revcomp(r.get("kmer", "")) in rec_kmers:
            checked += 1
            labs = r.get("selection_labels")
            if labs is None:
                raise Missing("selection_labels")
            if col not in labs:
                return FAIL, "mirror context %s not selected by %s" % (r["kmer"], col)
            ss = r.get("selection_strands")
            if ss is None:
                raise Missing("selection_strands")
            if ss[labs.index(col)] not in ("reverse_complement", "both"):
                return FAIL, "mirror context selected with strand %r" % ss[labs.index(col)]
    if checked < a.get("min", 1):
        # a withheld list or a cut can hide the mirrors, as in predicate_selected
        detail = "%d mirror context(s) returned, at least %d expected (the record holds the window's reverse " \
                 "complement in a longer inverted copy)" % (checked, a.get("min", 1))
        return _hidden_outcome(sc, ans, i, e, detail + "; the list is not complete (%s)", detail)
    return PASS, "%d mirror context(s) selected through their reverse complement" % checked


def check_motif(sc, ans, a, arm):
    e = _entry(ans, a.get("pattern", 0))
    m = e.get("motif")
    if m is None:
        raise Missing("motif")
    if "selected" in a and m.get("selected") != a["selected"]:
        return FAIL, "motif.selected %r, expected %r (decided_by %r)" % (m.get("selected"), a["selected"], m.get("decided_by"))
    if "decided_by_in" in a and m.get("decided_by") not in a["decided_by_in"]:
        return FAIL, "motif.decided_by %r, expected one of %r" % (m.get("decided_by"), a["decided_by_in"])
    if m.get("decided_by") == "tested_contexts" and m.get("labels_absent") is not None:
        return FAIL, "tested_contexts with labels_absent %r (must be null)" % m.get("labels_absent")
    if m.get("decided_by") == "no_instance" and m.get("labels_present") != []:
        return FAIL, "no_instance with labels_present %r" % m.get("labels_present")
    if "label_present" in a:
        cols = [x.get("column") for x in (m.get("labels_present") or [])]
        if a["label_present"] not in cols:
            return FAIL, "label %s not in labels_present %r" % (a["label_present"], cols)
    return PASS, "motif selected %r decided_by %r" % (m.get("selected"), m.get("decided_by"))


def check_clamped(sc, ans, a, arm):
    body = _body(ans)
    where = a.get("where", "limits")
    if where == "strategy":
        lst = _get(body, "strategy.clamped")
    else:
        lst = _get(body, "limits.clamped")
    for c in lst:
        if c.get("field") == a["field"]:
            if c.get("requested") != a["requested"] or c.get("effective") != a["effective"]:
                return FAIL, "clamp %r, expected requested %r effective %r" % (c, a["requested"], a["effective"])
            return PASS, "clamped %s %s -> %s" % (a["field"], a["requested"], a["effective"])
    return FAIL, "%s not in clamped %r" % (a["field"], lst)


# ---------------------------------------------------------------------------------------------- the time budget
def _number(v: Any) -> Optional[float]:
    return float(v) if isinstance(v, (int, float)) and not isinstance(v, bool) else None


# The pieces a walk may pass after its deadline before a reading of the clock stops it (SPEC-labeled-traversal-core
# §6.8; walker.cpp run(), run_level(), checkpoint(), pacing()). At a checked depth (>= 1) the deadline is read before
# every head and before every chunk of a paced read: at most one read (or chunk) and the head piece after it.
_CHECKED_PIECES = 2
# The depth-0 level of each requested arm is never checked (checkpoint: depth_ > 0) and its reads are unpaced
# (pacing(depth_ > 0 ? WALK : NONE)): its level fetch (under a request budget up to 3 calls for the root's at most
# 4 successor keys: 1, 2, then the rest), its lookahead (paced at depth 0 too: one chunk or chain piece past the
# deadline) and two head pieces (the root's successor enumeration before the fetch, its processing after it).
_DEPTH0_PIECES_PER_ARM = 6
# After the last reading that found no head below the radius: the radius level (its heads are only ended, with no
# reading) and the seed's finalisation.
_TAIL_PIECES = 2
# A 503 `deadline` of /pattern or /resolve comes at the deadline's next reading (SPEC-pattern-search §7.6,
# SPEC-labeled-traversal-core §4.5); the client's wall time also holds the network, the front and any wait for a
# request-pool thread, so a 503 is judged "long past" its budget only beyond max(1 s, the budget again).
_LATE_503_MIN_MS = 1000.0


def _budget_503(sc, ans, client, a=None):
    """budget_respected on an allowed 503 `deadline` of /pattern or /resolve: the client's wall time (the runner's
    `client`) against the EFFECTIVE budget, min(the requested budget, the route's cap the builder states in the
    check's args, cap_ms: /resolve resolve.time_budget.max_time_ms, /pattern caps.time_budget_ms of the capabilities
    snapshot; a host lowers a larger budget to its cap, SPEC-labeled-traversal-core §4.5, SPEC-pattern-search §4.5;
    A2-03 of the review of lane S5b). FAIL when it came before that budget, on both routes: the SPECs answer it
    only when the answer cannot be written by the budget (the host throws it once the elapsed time reaches the
    budget, Deadline::respond_expired; /resolve's finalisation reserve ends the work earlier, not the answer), and
    the deadline starts at the parsed body, inside the client's interval (SPEC-pattern-search §7.6,
    SPEC-labeled-traversal-core §4.5); FINDING when it came long past it (an unpolled piece of that length, which the
    SPECs' deadline rule answers at its next reading); else INCONCLUSIVE (nothing in the body to judge)."""
    asked = _stated_budget(sc)
    wall = _number((client or {}).get("wall_ms"))
    code = (ans.get("error") or {}).get("code")
    if asked is None or wall is None or code != "deadline":
        return INCONCLUSIVE, "not testable: http %s (%s), no client wall time against a budget" % (
            ans.get("http"), code)
    cap = _number((a or {}).get("cap_ms"))
    budget = min(asked, cap) if cap is not None else asked
    named = "%.0f ms budget" % budget if budget == asked else \
        "%.0f ms budget (the requested %.0f ms lowered to the host's cap)" % (budget, asked)
    late = wall - budget
    if wall < budget - 1.0:
        return FAIL, ("503 deadline at a client wall time of %.0f ms, before the %s: the answer is a 503 only when it "
                      "cannot be written by the budget, and the client's interval holds the host's" % (wall, named))
    if late > max(_LATE_503_MIN_MS, budget):
        return FINDING, ("503 deadline delivered %.0f ms after the %s (client wall %.0f ms): the deadline is read "
                         "where the work and the delivery poll it, so an unpolled piece of about that length ran (the "
                         "wall also holds the network and any wait for a request thread)" % (late, named, wall))
    return INCONCLUSIVE, "503 deadline at a client wall time of %.0f ms for a %s (within max(1 s, the budget) of " \
                         "it): nothing in the body to judge" % (wall, named)


def _seed_paths_max(sc, r, i):
    """The longest path of the requested arms of seed result |r| (bp; 0 without paths)."""
    best = 0
    for side in _arms_of(sc, i):
        for p in ((r.get("arms") or {}).get(side) or {}).get("paths") or []:
            if isinstance(p.get("length_bp"), int):
                best = max(best, p["length_bp"])
    return best


def _traverse_time_stop(r: dict, budget: float) -> Tuple[Optional[float], Optional[float]]:
    """(when the walk saw its time stop at the latest, when it saw it) on the seed's clock: the former the smallest
    of timing.elapsed_ms, resource_stop.used and budget + stop_after_deadline_ms (walker.cpp note_stop clamps that
    at 0, so a premature stop shows 0 there and its time in `used` or the elapsed time), the latter budget +
    stop_after_deadline_ms when positive, else resource_stop.used (None when unknown)."""
    t = r.get("timing") or {}
    dl = t.get("deadline") or {}
    rs = r.get("resource_stop") or {}
    el = _number(t.get("elapsed_ms"))
    sad = _number(dl.get("stop_after_deadline_ms"))
    used = _number(rs.get("used")) if rs.get("resource") == "time" else None
    known = [x for x in (el, used, (budget + sad) if sad is not None else None) if x is not None]
    stop = budget + sad if sad is not None and sad > 0 else used
    return (min(known) if known else None), stop


def check_budget_respected(sc, ans, a, arm, client=None):
    """The server's time against the request's effective time budget, each route by its SPEC's rule and the
    answer's own deadline record; no fixed allowance (the capabilities state no bound on one uninterruptible
    piece: deadline_check.max_uninterruptible_ms is null). /pattern and /resolve: a 200 whose timing.elapsed_ms
    exceeds limits.time_budget_ms broke the SPEC (past the budget they answer 503); an allowed 503 `deadline` long
    past the budget by the client's wall time (_budget_503) is a finding too. /traverse (SPEC-labeled-traversal-core
    §6.8, walker.cpp), per walked seed: the deadline is the budget on the seed's clock (time_exceeded: elapsed >=
    budget), never read at depth 0, the depth-0 reads unpaced, so
      - a time stop before the budget is a finding (the walker stops only once the budget has passed);
      - a budget that passed in the seed phase stops the walk at the first depth-1 checkpoint: a path longer than
        1 bp is a finding;
      - a time stop at a checked depth (cap_trigger.at_bp >= 2) comes within one read and one head piece after
        max(budget, seed phase), each at most the seed's longest piece (timing.deadline.longest_piece, and at most
        the time since the seed phase); a stop at depth 1 (at_bp <= 1) within the depth-0 levels' pieces of the
        requested arms (_DEPTH0_PIECES_PER_ARM each); later is a finding;
      - a walk with no time stop past max(budget, seed phase) by more than the pieces that may follow the deadline
        before a reading that would stop a head below the radius (_CHECKED_PIECES, or the depth-0 pieces when no
        path passed 1 bp, plus _TAIL_PIECES: the radius level and the finalisation) is a finding.
    A breach is a FINDING (the request still answered); a stop late within its pieces, or a seed phase longer than
    the budget, is a performance fact stated in the detail (the report's timing section)."""
    body = _body(ans)
    route = sc.get("route", "")
    if route == "POST /traverse":
        budget = _effective_budget(sc, body)
        if budget is None:
            raise Missing("strategy.bounds.time_budget_ms")
        res = _get(body, "results")
        over, facts, seen = [], [], 0
        for i in range(len(res)):
            r = _result(ans, i)
            if "arms" not in r:
                continue   # no walk: the seed phase alone, which the budget does not govern
            t = r.get("timing")
            if not isinstance(t, dict) or _number(t.get("elapsed_ms")) is None:
                raise Missing("results[%d].timing.elapsed_ms" % i)
            seen += 1
            el = _number(t["elapsed_ms"])
            sp = _number(t.get("seed_phase_ms"))
            if sp is None:
                raise Missing("results[%d].timing.seed_phase_ms" % i)
            dl = t.get("deadline") or {}
            lp = dl.get("longest_piece") or {}
            lp_ms = _number(lp.get("ms")) or 0.0
            # every piece after the deadline lies after the seed phase, so none is longer than the time since then
            piece_ms = min(lp_ms, max(0.0, el - sp))
            piece = "%s of %.0f ms (%s rows%s)" % (lp.get("kind"), lp_ms, lp.get("rows"),
                                                   ", %.0f ms since the seed phase" % (el - sp) if piece_ms < lp_ms else "")
            arms = _arms_of(sc, i)
            deadline = max(budget, sp)
            longest = _seed_paths_max(sc, r, i)
            if sp >= budget:
                facts.append("seed %d: the seed phase (%.0f ms) outlasted the budget, which does not govern it" % (i, sp))
                if longest > 1:
                    over.append("seed %d: the budget passed in the seed phase (%.0f ms of %.0f), so the first depth-1 "
                                "checkpoint stops the walk, yet a path runs %d bp (at most 1 expected)" % (
                                    i, sp, budget, longest))
            time_stop = dl.get("stopped_by") == "time_budget" or (r.get("resource_stop") or {}).get("resource") == "time"
            if time_stop:
                latest, stop = _traverse_time_stop(r, budget)
                if latest is not None and latest < budget - 1.0:
                    over.append("seed %d: a time stop by %.0f ms, before its %.0f ms budget (the walker stops only once "
                                "the budget has passed)" % (i, latest, budget))
                    continue
                ats = [((r.get("arms") or {}).get(x) or {}).get("cap_trigger") or {} for x in arms]
                ats = [c.get("at_bp") for c in ats if c.get("reason") == "time_budget" and isinstance(c.get("at_bp"), int)]
                at = min(ats) if ats else None
                if stop is None:
                    raise Missing("results[%d].timing.deadline.stop_after_deadline_ms" % i)
                late = stop - deadline
                if at is not None and at >= 2:
                    allow, where = _CHECKED_PIECES * piece_ms, "at a checked depth (%d bp)" % at
                else:
                    allow = _DEPTH0_PIECES_PER_ARM * len(arms) * piece_ms
                    where = "at depth 1 (after the unpaced, unchecked depth-0 level of %d arm(s))" % len(arms)
                if late > allow:
                    over.append("seed %d: stopped %s %.0f ms after max(budget %.0f, seed phase %.0f), beyond the %.0f ms "
                                "its pieces allow (longest piece a %s)" % (i, where, late, budget, sp, allow, piece))
                elif late > 0:
                    facts.append("seed %d: stopped %s %.0f ms after the deadline, within its pieces (longest a %s; "
                                 "elapsed %.0f ms)" % (i, where, late, piece, el))
            elif el > deadline:
                walked_on = longest >= 2   # the depth-1 reading came before the deadline
                pieces = (_CHECKED_PIECES if walked_on else _DEPTH0_PIECES_PER_ARM * len(arms)) + _TAIL_PIECES
                allow = pieces * piece_ms
                complete = all(((r.get("arms") or {}).get(x) or {}).get("status") == "complete" for x in arms)
                if el - deadline > allow:
                    over.append("seed %d: elapsed %.0f ms, %.0f ms past max(budget %.0f, seed phase %.0f) without a time "
                                "stop (%s), beyond the %d pieces (%.0f ms; longest a %s) that may follow the deadline "
                                "before a reading that stops a head below the radius" % (
                                    i, el, el - deadline, budget, sp, "every arm complete" if complete else
                                    "arms ended otherwise", pieces, allow, piece))
                else:
                    facts.append("seed %d: elapsed %.0f ms past the budget without a time stop, within %d pieces (%s)" % (
                        i, el, pieces, "a level at the radius is never stopped" if complete else "arms ended otherwise"))
        if over:
            return FINDING, "; ".join(over[:4])
        if not seen:
            return INCONCLUSIVE, "no walked seed: nothing the budget governs"
        return PASS, "%d walked seed(s) kept the budget %.0f ms by the deadline record%s" % (
            seen, budget, ("; performance: " + "; ".join(facts)) if facts else "")
    if ans.get("http") == 503:
        return _budget_503(sc, ans, client, a)
    budget = _effective_budget(sc, body)
    if budget is None:
        raise Missing("limits.time_budget_ms")
    el = _number(_get(body, "timing.elapsed_ms"))
    if el is None:
        raise Missing("timing.elapsed_ms")
    if el > budget:
        stops = []
        for i, e in enumerate(body.get("patterns") or []):
            if isinstance(e, dict) and e.get("stop"):
                stops.append("p%d stop %s/%s" % (i, e["stop"].get("phase"), e["stop"].get("reason")))
        if body.get("stop"):
            stops.append("stop %s/%s" % (body["stop"].get("phase"), body["stop"].get("reason")))
        return FINDING, "a 200 after %.0f ms, past the time budget %.0f ms (%s)" % (
            el, budget, ", ".join(stops) or "no stop stated")
    return PASS, "elapsed %.0f ms within the budget %.0f ms" % (el, budget)


# ---------------------------------------------------------------------------------------------- pattern invariants
def _sum_counts(parts: List[dict]) -> Tuple[Optional[int], str]:
    """§7.4's algebra over counts (without lower/upper detail: value and relation)."""
    rels = [p["relation"] for p in parts]
    if not parts:
        return 0, "exact"
    if all(r == "unknown" for r in rels):
        return None, "unknown"
    if any(r in ("unknown", "at_least") for r in rels):
        return sum((p.get("lower", p.get("value")) or 0) if p["relation"] == "bounds" else (p.get("value") or 0)
                   for p in parts), "at_least"
    if any(r == "bounds" for r in rels):
        return sum(p.get("lower", p.get("value")) or 0 for p in parts), "bounds"
    return sum(p["value"] for p in parts), "exact"


def _check_count_shape(c: dict, where: str) -> Optional[str]:
    rel = c.get("relation")
    if rel not in ("exact", "at_least", "bounds", "unknown"):
        return "%s: relation %r" % (where, rel)
    if (c.get("value") is None) != (rel == "unknown"):
        return "%s: value %r with relation %s" % (where, c.get("value"), rel)
    if rel == "bounds" and not (c.get("lower", -1) <= c.get("upper", -2)):
        return "%s: bounds without lower <= upper" % where
    return None


def _pattern_entry_problems(sc: dict, i: int, e: dict, limits: dict, sticky: Optional[dict]) -> List[str]:
    probs: List[str] = []
    req = sc["request"]
    try:
        kind, pat = _request_pattern(sc, i)
    except (Missing, KeyError, IndexError, TypeError):
        return probs
    if "error" in e:
        if "counts" in e or "results" in e:
            probs.append("p%d: refused slot with counts/results" % i)
        return probs
    if "counts" not in e:
        return probs
    table = req.get("genetic_code", 1) if kind == "protein" else 1
    L = _pattern_length(kind, pat)
    mode = e.get("mode", req.get("mode", "all_or_count"))
    counts = e["counts"]
    for name, c in counts.items():
        if isinstance(c, dict) and "relation" in c:
            m = _check_count_shape(c, "p%d counts.%s" % (i, name))
            if m:
                probs.append(m)
    m_max = req.get("max_mismatches", 0) if isinstance(req.get("max_mismatches", 0), int) else 0
    if L <= K and "contexts" in counts:
        cc = counts["contexts"]
        bo = cc.get("by_offset", {})
        offsets = [str(K - L)] if e.get("scope") == "suffix" else [str(p) for p in range(0, K - L + 1)]
        if sorted(bo.keys(), key=int) != sorted(offsets, key=int):
            probs.append("p%d: by_offset keys %s, expected %s" % (i, sorted(bo.keys(), key=int)[:5], offsets[:5]))
        for part in ("by_offset", "by_strand", "by_orientation"):
            if part in cc and cc.get("relation") == "exact":
                v, r = _sum_counts(list(cc[part].values()))
                if r == "exact" and v != cc["value"]:
                    probs.append("p%d: %s sums to %s, total %s" % (i, part, v, cc["value"]))
        if m_max > 0:
            bm = cc.get("by_mismatches")
            if bm is None:
                probs.append("p%d: by_mismatches absent with max_mismatches %d" % (i, m_max))
            else:
                if sorted(bm.keys(), key=int) != [str(t) for t in range(m_max + 1)]:
                    probs.append("p%d: by_mismatches keys %s" % (i, sorted(bm.keys())))
                v, r = _sum_counts([bm[str(t)] for t in range(m_max + 1) if str(t) in bm])
                if r != cc["relation"] and not (r == "at_least" and cc["relation"] == "at_least"):
                    probs.append("p%d: levels sum to relation %s, total %s" % (i, r, cc["relation"]))
                if r == "exact" and v != cc.get("value"):
                    probs.append("p%d: levels sum to %s, total %s" % (i, v, cc.get("value")))
                stop = e.get("stop")
                if stop and "mismatch_level" in stop:
                    t0 = stop["mismatch_level"]
                    for t in range(m_max + 1):
                        rel = bm.get(str(t), {}).get("relation")
                        if t < t0 and rel == "unknown":
                            probs.append("p%d: level %d unknown below the stop's level %d" % (i, t, t0))
                        if t > t0 and rel != "unknown":
                            probs.append("p%d: level %d %s above the stop's level %d" % (i, t, rel, t0))
                sbm = (e.get("work") or {}).get("steps_by_mismatches")
                if sbm is not None and sum(sbm.values()) != e["work"].get("steps"):
                    probs.append("p%d: steps_by_mismatches sum %d != steps %d" % (i, sum(sbm.values()), e["work"].get("steps")))
    # results
    if mode == "count":
        for f in ("results", "withheld", "returned", "cut"):
            if f in e:
                probs.append("p%d: mode count with %s" % (i, f))
        if e.get("retrieval_complete") is not False:
            probs.append("p%d: mode count with retrieval_complete %r" % (i, e.get("retrieval_complete")))
        return probs
    res = e.get("results")
    if res is None:
        probs.append("p%d: retrieval mode without results" % i)
        return probs
    if e.get("returned") != len(res):
        probs.append("p%d: returned %s != len(results) %d" % (i, e.get("returned"), len(res)))
    if e.get("withheld") and res:
        probs.append("p%d: withheld %r with %d results" % (i, e["withheld"], len(res)))
    if mode != "partial" and e.get("cut"):
        probs.append("p%d: cut %r outside partial" % (i, e["cut"]))
    keys = []
    for j, r in enumerate(res):
        strand = r.get("strand", r.get("orientation"))
        if "kmer" in r:
            off = r.get("offset", -1)
            if not (0 <= off <= K - L) or r.get("instance") != r["kmer"][off:off + L]:
                probs.append("p%d r%d: instance %r is not kmer[%s:%s]" % (i, j, r.get("instance"), off, off + L))
                continue
            inst = r["instance"]
            lev = r.get("mismatches", 0)
            if m_max > 0 or "substitutions" in r:
                subs = expected_substitutions(kind, pat, inst, strand, table) if strand != "=" else []
                if subs is None:
                    probs.append("p%d r%d: instance %s is no instance at any level (SPEC-pattern-search §21.3: a stop "
                                 "codon never substitutes a residue, and X and * hold no substitution)" % (i, j, inst))
                elif strand != "=":
                    if len(subs) != lev:
                        probs.append("p%d r%d: level %d but %d substitutions recomputed" % (i, j, lev, len(subs)))
                    elif r.get("substitutions") != subs:
                        probs.append("p%d r%d: substitutions %r, recomputed %r" % (i, j, r.get("substitutions"), subs))
                if lev > m_max:
                    probs.append("p%d r%d: level %d above max_mismatches %d" % (i, j, lev, m_max))
            elif strand in ("+", "-", "forward", "reverse") and not _instance_matches(kind, pat, inst, strand, table):
                probs.append("p%d r%d: instance %s does not match the pattern on %s" % (i, j, inst, strand))
            keys.append((lev, r.get("node", 0), off, _ORIENT_RANK.get(strand, 3)))
        elif "sequence" in r:
            if r.get("instance") not in (None, r["sequence"]) or len(r["sequence"]) != L:
                probs.append("p%d r%d: path sequence of %d bases, instance %r" % (i, j, len(r["sequence"]), r.get("instance", "")[:10]))
            elif strand in ("+", "-") and not _instance_matches(kind, pat, r["sequence"], strand, table):
                probs.append("p%d r%d: path %s... does not spell the pattern on %s" % (i, j, r["sequence"][:20], strand))
            if "nodes" in r and len(r["nodes"]) != L - K + 1:
                probs.append("p%d r%d: %d nodes for L=%d" % (i, j, len(r["nodes"]), L))
            keys.append((0, (r.get("nodes") or [0])[0], _ORIENT_RANK.get(strand, 3), r["sequence"]))
    if keys != sorted(keys):
        probs.append("p%d: results not in answer order" % i)
    # a context is one k-mer at an offset in an orientation, listed at one level (SPEC-pattern-search §7.1,
    # §21.4), a path its nodes in an orientation (§12.1): each is listed once, whatever the order says of equal keys
    # (a copy in place of a distinct result keeps the list's length, which counts_consistent compares)
    seen: Dict[tuple, int] = {}
    for j, r in enumerate(res):
        o = r.get("strand", r.get("orientation"))
        if "kmer" in r:
            ident = ("context", r.get("node"), r.get("kmer"), r.get("offset"), o)
        elif "sequence" in r:
            ident = ("path", tuple(r.get("nodes") or ()), r.get("sequence"), o)
        else:
            continue
        if ident in seen:
            probs.append("p%d r%d: the result of r%d listed again (%s %s)" % (
                i, j, seen[ident], o, ("offset %s of %s" % (r.get("offset"), r.get("kmer"))) if "kmer" in r
                else "path %s..." % str(r.get("sequence"))[:20]))
        else:
            seen[ident] = j
    if e.get("retrieval_complete"):
        name = _primary_count_name(e, req)
        c = counts.get(name)
        if not c or c.get("relation") != "exact" or c.get("value") != len(res):
            probs.append("p%d: retrieval_complete with %s %r and %d results" % (i, name, c and (c.get("relation"), c.get("value")), len(res)))
    stop = e.get("stop")
    if _note_stop(sc, e) == "budget" and "low_complexity_pattern" in (e.get("notes") or []):
        # the engine diagnoses the note only after a search that completed or stopped at a threshold
        # (SPEC-pattern-search §7.8); a later phase's stop leaves it as it is
        probs.append("p%d: low_complexity_pattern beside a budget stop %s/%s of the search" % (
            i, stop.get("phase"), stop.get("reason")))
    # SPEC-pattern-search §7.5: a stop of the graph search shortens a list only with a cut (partial) or a withheld
    # list (all_or_count)
    if isinstance(stop, dict) and stop.get("phase") in ("discovery", "mask_scan", "extraction", "extension", "selection") and \
            not e.get("withheld") and not e.get("cut"):
        probs.append("p%d: a %s/%s stop with neither a cut nor a withheld list" % (i, stop.get("phase"), stop.get("reason")))
    return probs


def _counts_consistent_problems(sc: dict, i: int, e: dict, limits: dict) -> List[str]:
    probs: List[str] = []
    req = sc["request"]
    if "counts" not in e or "error" in e:
        return probs
    mode = e.get("mode", req.get("mode", "all_or_count"))
    if mode == "count":
        return probs
    res = e.get("results") or []
    name = _primary_count_name(e, req)
    c = e["counts"].get(name)
    if not c:
        return probs
    value, rel = c.get("value"), c.get("relation")
    is_long = "contexts" not in e["counts"]
    cap = limits.get("max_paths") if (is_long and name in ("paths", "supported_paths", "selected")) else limits.get("max_contexts")
    if cap is None:
        return probs
    stop = e.get("stop")
    if len(res) > cap:
        probs.append("p%d: %d results above the threshold %d" % (i, len(res), cap))
    if mode == "all_or_count":
        if e.get("withheld") is None and rel == "exact" and value != len(res) and not stop:
            probs.append("p%d: all_or_count released %d of exact %d" % (i, len(res), value))
        if rel == "exact" and value is not None and value > cap and e.get("withheld") is None:
            probs.append("p%d: all_or_count exact %d above %d not withheld" % (i, value, cap))
    if mode == "partial" and rel == "exact" and not stop and e.get("withheld") is None:
        if value > cap:
            want = "max_paths" if (is_long and cap == limits.get("max_paths")) else "max_contexts"
            if (e.get("cut") or {}).get("reason") not in (want, "max_memory") or (len(res) != cap and (e.get("cut") or {}).get("reason") != "max_memory"):
                probs.append("p%d: partial exact %d above %d: cut %r, %d returned" % (i, value, cap, e.get("cut"), len(res)))
        elif e.get("cut") is None and len(res) != value:
            probs.append("p%d: partial exact %d below the cap but %d returned, no cut" % (i, value, len(res)))
    expected_rc = (rel == "exact" and value == len(res) and e.get("withheld") is None)
    if e.get("retrieval_complete") and not expected_rc:
        probs.append("p%d: retrieval_complete true without an exact count equal to the list" % i)
    return probs


# the counts a sticky stop leaves as they are (no population count: §19.7 `tested`, the decisions a selection pass
# made; `estimate`, no bound), as compare.NOT_POPULATION
_STICKY_NOT_POPULATION = ("tested", "estimate")


def _known_counts(counts: Any, path: str = "") -> List[str]:
    """Every count object of an entry's |counts| (each part: suffix, by_offset, by_strand / by_orientation,
    by_mismatches) whose relation is not unknown, as "path relation value"."""
    out: List[str] = []
    if isinstance(counts, dict):
        if "relation" in counts and counts.get("relation") != "unknown":
            out.append("%s %s %s" % (path or "count", counts.get("relation"), counts.get("value")))
        for k, v in counts.items():
            if str(k) not in _STICKY_NOT_POPULATION and isinstance(v, dict):
                out += _known_counts(v, (path + "." if path else "") + str(k))
    return out


def _sticky_problems(sc: dict, i: int, e: dict, first: dict) -> List[str]:
    """What entry i, a later pattern after |first|, an earlier pattern's budget stop of the search, states otherwise
    than SPEC-pattern-search §7.6 has it (one budget per request, the stop sticky): the same stop in phase discovery
    (with max_mismatches above 0 at level 0, §21.7), every count unknown (by_mismatches included, §21.6), work zero
    (steps_by_mismatches included); as the oracle's sticky_stop / stop.sticky_unknown (SD-013)."""
    req = sc.get("request") if isinstance(sc.get("request"), dict) else {}
    m = req.get("max_mismatches")
    want = {"phase": "discovery", "reason": first.get("reason")}
    if isinstance(m, int) and not isinstance(m, bool) and m > 0:
        want["mismatch_level"] = 0
    s = e.get("stop")
    out = []
    if s != want:
        out.append("p%d: after a %s stop of an earlier pattern, stop %r (the sticky stop is %r: "
                   "SPEC-pattern-search §7.6)" % (i, first.get("reason"), s, want))
    known = _known_counts(e.get("counts"))
    if known:
        out.append("p%d: after a %s stop of an earlier pattern, %d count(s) not unknown (%s; every count is unknown: "
                   "§7.6)" % (i, first.get("reason"), len(known), known[0]))
    work = e.get("work") if isinstance(e.get("work"), dict) else {}
    busy = ["%s %s" % (k, v) for k, v in work.items() if isinstance(v, (int, float)) and not isinstance(v, bool) and v]
    busy += ["steps_by_mismatches.%s %s" % (k, v) for k, v in (work.get("steps_by_mismatches") or {}).items()
             if isinstance(v, (int, float)) and v]
    if busy:
        out.append("p%d: after a %s stop of an earlier pattern, work %s (zero: §7.6)" % (
            i, first.get("reason"), ", ".join(busy[:3])))
    return out


def check_counts_consistent(sc, ans, a, arm):
    body = _body(ans)
    pats = _get(body, "patterns")
    limits = body.get("limits") or {}
    which = [a["pattern"]] if "pattern" in a else range(len(pats))
    probs: List[str] = []
    for i in which:
        probs += _counts_consistent_problems(sc, i, pats[i], limits)
    # one budget per request: a budget stop is sticky for the later patterns (SPEC §7.6; _sticky_problems); a refused
    # slot and a peptide answered without a search (note no_stop_codon, §12.2) are no later search
    first = None
    for i, e in enumerate(pats):
        s = e.get("stop")
        if first is not None and "counts" in e and "no_stop_codon" not in (e.get("notes") or []):
            probs += _sticky_problems(sc, i, e, first)
        if first is None and s and s.get("reason") in _BUDGET_STOP_REASONS and s.get("phase") in (
                "discovery", "mask_scan", "extraction", "extension"):
            first = s
    return (FAIL, "; ".join(probs[:6])) if probs else (PASS, "modes, lists and relations consistent")


def _pattern_invariants(sc, ans, a, arm):
    body = _body(ans)
    if ans.get("http") != 200:
        return PASS, "no answer to check (http %s)" % ans.get("http")
    pats = _get(body, "patterns")
    if len(pats) != len(sc["request"].get("patterns", [])):
        return FAIL, "%d entries for %d patterns" % (len(pats), len(sc["request"]["patterns"]))
    probs: List[str] = []
    undecided: List[str] = []
    limits = body.get("limits") or {}
    for i, e in enumerate(pats):
        req_id = sc["request"]["patterns"][i].get("id") if isinstance(sc["request"]["patterns"][i], dict) else None
        if e.get("id") != req_id:
            probs.append("p%d: id %r, request %r" % (i, e.get("id"), req_id))
        probs += _pattern_entry_problems(sc, i, e, limits, None)
        try:
            probs += _pattern_reason_problems(sc, body, i, e)
        except Undecidable as u:
            undecided.append(str(u))
    if probs:
        return FAIL, "; ".join(probs[:6])
    if undecided:
        raise Undecidable(undecided[0])
    return PASS, "%d entries consistent" % len(pats)


# ---------------------------------------------------------------------------------------------- resolve checks
def _resolve_label(body: dict, name: str) -> dict:
    labs = [l for l in _get(body, "labels") if l.get("label") == name]
    if not labs:
        raise Missing("label %s" % name)
    return labs[0]


def check_resolve_runs(sc, ans, a, arm, client=None):
    body = _body(ans)
    lab = _resolve_label(body, a["label"])
    n = len([l for l in body["labels"] if l.get("label") == a["label"]])
    if n != 1:
        return FAIL, "label %s listed %d times" % (a["label"], n)
    runs = lab.get("runs")
    if not isinstance(runs, list):
        raise Missing("runs")
    num = body.get("num_kmers")
    ks = lab.get("kmers_supported")
    stop = body.get("stop")
    if stop:
        bad = _resolve_stop_problem(sc, body, client)
        if bad:
            return FAIL, "label %s: %s" % (a["label"], bad)
    if (a.get("prefix_if_stopped") or a.get("unless_stopped")) and stop:
        # a stop answers the resolve of the prefix of resolved_kmers k-mers (SPEC §4.5): the record's own label
        # supports every k-mer of it; a stop before the first k-mer shows nothing
        x = stop.get("resolved_kmers")
        if not x:
            return INCONCLUSIVE, "stopped before the first k-mer (resolved_kmers %r): no prefix to check" % x
        if a.get("full") and (ks != x or runs != [[0, x]] or body.get("graph_runs") != [[0, x]]):
            return FAIL, "stopped at %s k-mers: label %s supports %s, runs %s, graph_runs %s (expected [0, %s))" % (
                x, a["label"], ks, runs[:4], (body.get("graph_runs") or [])[:4], x)
        return PASS, "stopped after %s of %s k-mers: the prefix fully supported by %s" % (
            x, stop.get("query_kmers"), a["label"])
    if stop:
        # the answer is the resolve of a prefix (SPEC-labeled-traversal-core §4.5), not of the query the expectation
        # is about: each argument is decided only where the prefix decides it
        return _resolve_prefix_profile(a, lab, runs, stop)
    outcome, detail = _resolve_profile(a, lab, runs, num, ks)
    if a.get("prefix_if_stopped") and not stop and outcome != FAIL:
        # the scenario's property is the stated stop (SPEC §4.5: a stop answers the resolve of a prefix): a complete
        # answer holds the profile but never exercised it
        return INCONCLUSIVE, "no stop happened: the prefix property was not tested (%s)" % detail
    return outcome, detail


def _resolve_prefix_profile(a: dict, lab: dict, runs: list, stop: dict) -> Tuple[str, str]:
    """resolve_runs' assertions about the whole query on a stopped answer, which resolves the prefix of the first
    resolved_kmers = x k-mers (SPEC-labeled-traversal-core §4.5: the profile of the prefix, a run that reaches the
    prefix's end cut there). An argument fails where the prefix contradicts it (a hole in a full support, a run end
    before the range, more support than allowed, too little left to reach a count), passes where the prefix alone
    proves it (a covered interval inside the prefix, a run end inside it, a minimum already reached) and is
    inconclusive otherwise; a stop before the first k-mer shows nothing."""
    label = a["label"]
    x = _number(stop.get("resolved_kmers"))
    if not x:
        return INCONCLUSIVE, "stopped before the first k-mer (resolved_kmers %r): nothing of label %s to check" % (
            stop.get("resolved_kmers"), label)
    x = int(x)
    qk = _number(stop.get("query_kmers"))
    rest = int(qk) - x if qk is not None and qk >= x else None
    ks = _number(lab.get("kmers_supported"))
    clean = [r for r in runs if isinstance(r, list) and len(r) == 2]
    head = "stopped at %d of %s k-mers: label %s" % (x, stop.get("query_kmers"), label)
    results: List[Tuple[str, str]] = []
    if a.get("full"):
        if ks != x or clean != [[0, x]]:
            results.append((FAIL, "%s supports %s, runs %s: the full support expected is broken in the prefix" % (
                head, lab.get("kmers_supported"), clean[:4])))
        else:
            results.append((INCONCLUSIVE, "%s supports the whole prefix; full support of the query untested" % head))
    if "covers" in a:
        lo, hi = a["covers"]
        if hi <= x:
            ok = any(r[0] <= lo and hi <= r[1] for r in clean)
            results.append((PASS if ok else FAIL, "%s: [%d, %d) %scovered in the prefix (runs %s)" % (
                head, lo, hi, "" if ok else "not ", clean[:6])))
        elif lo < x and not any(r[0] <= lo and r[1] == x for r in clean):
            results.append((FAIL, "%s: no run covers [%d, %d) of the prefix, so none covers [%d, %d) (runs %s)" % (
                head, lo, x, lo, hi, clean[:6])))
        else:
            results.append((INCONCLUSIVE, "%s: [%d, %d) reaches past the prefix" % (head, lo, hi)))
    if "run_end_in" in a:
        lo, hi = a["run_end_in"]
        first = next((r for r in clean if r[0] == 0), None)
        if first is None:
            results.append((FAIL, "%s: no run starts at the first k-mer (runs %s)" % (head, clean[:6])))
        elif first[1] < x:
            ok = lo <= first[1] <= hi
            results.append((PASS if ok else FAIL, "%s: the run [0, %d) ends inside the prefix, %s [%d, %d]" % (
                head, first[1], "in" if ok else "not in", lo, hi)))
        elif x > hi:
            results.append((FAIL, "%s: the run from 0 reaches the prefix's end %d, past [%d, %d]" % (head, x, lo, hi)))
        else:
            results.append((INCONCLUSIVE, "%s: the run from 0 reaches the prefix's end %d, before [%d, %d] ends" % (
                head, x, lo, hi)))
    if "run_start_in" in a:
        lo, hi = a["run_start_in"]
        if x > hi and not any(r[1] == x and lo <= r[0] <= hi for r in clean):
            results.append((FAIL, "%s: no run [b, %d) with b in [%d, %d], so none reaches the query's end from "
                                  "there (runs %s)" % (head, x, lo, hi, clean[:6])))
        else:
            results.append((INCONCLUSIVE, "%s: a run to the query's end is past the prefix" % head))
    if "runs" in a:
        want = [[r[0], min(r[1], x)] for r in a["runs"] if r[0] < x]
        if clean != want:
            results.append((FAIL, "%s: runs %s, the expected %s cut at the prefix are %s" % (
                head, clean[:6], a["runs"][:6], want[:6])))
        else:
            results.append((INCONCLUSIVE, "%s: runs %s as expected in the prefix" % (head, clean[:6])))
    for arg, bad, good in (("kmers", lambda v: ks > v or (rest is not None and ks + rest < v), None),
                           ("max_kmers", lambda v: ks > v, None),
                           ("min_kmers", lambda v: rest is not None and ks + rest < v, lambda v: ks >= v)):
        if arg not in a or ks is None:
            continue
        if bad(a[arg]):
            results.append((FAIL, "%s supports %s k-mers of the prefix, %s k-mers left: %s %s cannot hold" % (
                head, lab.get("kmers_supported"), rest, arg, a[arg])))
        elif good is not None and good(a[arg]):
            results.append((PASS, "%s supports %s k-mers of the prefix, %s %s reached" % (
                head, lab.get("kmers_supported"), arg, a[arg])))
        else:
            results.append((INCONCLUSIVE, "%s supports %s k-mers of the prefix; %s %s untested" % (
                head, lab.get("kmers_supported"), arg, a[arg])))
    for want in (FAIL, INCONCLUSIVE):
        hit = [d for o, d in results if o == want]
        if hit:
            return want, "; ".join(hit[:3])
    if results:
        return PASS, "; ".join(d for _, d in results[:3])
    return INCONCLUSIVE, "%s: nothing of the query's profile asked" % head


def _resolve_profile(a: dict, lab: dict, runs: list, num: Any, ks: Any) -> Tuple[str, str]:
    """resolve_runs' assertions on a complete profile of label a["label"]: full support, covered intervals, run
    ends and starts, exact runs, the supported k-mers."""
    if a.get("full"):
        if ks != num or runs != [[0, num]]:
            return FAIL, "label %s: %s of %s k-mers, runs %s (expected one run [0, %s])" % (a["label"], ks, num, runs[:4], num)
    for (lo, hi) in [a["covers"]] if "covers" in a else []:
        if not any(r[0] <= lo and hi <= r[1] for r in runs):
            return FAIL, "label %s: no run covers [%d, %d) (runs %s)" % (a["label"], lo, hi, runs[:6])
    if "run_end_in" in a:
        lo, hi = a["run_end_in"]
        if not any(r[0] == 0 and lo <= r[1] <= hi for r in runs):
            return FAIL, "label %s: no run [0, e) with e in [%d, %d] (runs %s)" % (a["label"], lo, hi, runs[:6])
    if "run_start_in" in a:
        lo, hi = a["run_start_in"]
        if not any(r[1] == num and lo <= r[0] <= hi for r in runs):
            return FAIL, "label %s: no run [b, %s) with b in [%d, %d] (runs %s)" % (a["label"], num, lo, hi, runs[:6])
    if "runs" in a and runs != a["runs"]:
        return FAIL, "label %s: runs %s, expected exactly %s" % (a["label"], runs[:6], a["runs"])
    if "kmers" in a and ks != a["kmers"]:
        return FAIL, "label %s supports %s k-mers, expected %s" % (a["label"], ks, a["kmers"])
    if "max_kmers" in a and ks > a["max_kmers"]:
        return FAIL, "label %s supports %s k-mers, expected at most %s" % (a["label"], ks, a["max_kmers"])
    if "min_kmers" in a and ks < a["min_kmers"]:
        return FAIL, "label %s supports %s k-mers, expected at least %s" % (a["label"], ks, a["min_kmers"])
    return PASS, "label %s: %s of %s k-mers, %d run(s)" % (a["label"], ks, num, len(runs))


def check_label_not_full(sc, ans, a, arm):
    body = _body(ans)
    lab = _resolve_label(body, a["label"])
    if lab.get("kmers_supported") == body.get("num_kmers") and body.get("num_kmers"):
        return FAIL, "unrelated label %s supports every k-mer" % a["label"]
    return PASS, "label %s supports %s of %s k-mers" % (a["label"], lab.get("kmers_supported"), body.get("num_kmers"))


def check_resolve_absent(sc, ans, a, arm):
    body = _body(ans)
    gr = _get(body, "graph_runs")
    if gr:
        return FAIL, "graph_runs %s, expected none" % gr[:4]
    for lab in body.get("labels") or []:
        if lab.get("kmers_supported"):
            return FAIL, "label %s supports %s k-mers of an absent sequence" % (lab.get("label"), lab["kmers_supported"])
    return PASS, "no k-mer in the graph"


def check_resolve_graph_runs(sc, ans, a, arm):
    body = _body(ans)
    gr = _get(body, "graph_runs")
    if gr != a["runs"]:
        return FAIL, "graph_runs %s, expected %s" % (gr[:6], a["runs"])
    return PASS, "graph_runs %s" % gr


# a discovery's max_labels where the request states none (traverse.cpp parses discover.max_labels with default 1000
# and no server cap; SPEC-labeled-traversal-core §4.1's example states the same)
_DISCOVER_MAX_LABELS = 1000


def _discover_max_labels(req: dict) -> Optional[float]:
    """The max_labels a discovery request runs with: the one it states, else the default."""
    d = req.get("discover") if isinstance(req, dict) else None
    if not isinstance(d, dict):
        return None
    return _number(d["max_labels"]) if "max_labels" in d else float(_DISCOVER_MAX_LABELS)


def _in_graph_kmers(body: dict) -> int:
    """The query's k-mers in the graph (the resolved prefix's under a stop): the lengths of graph_runs summed."""
    return sum(r[1] - r[0] for r in body.get("graph_runs") or [] if isinstance(r, list) and len(r) == 2)


def check_resolve_discover_contains(sc, ans, a, arm):
    """The discovery lists the scenario's source label (SPEC-labeled-traversal-core §4.2: the labels of the query's
    in-graph k-mers; more than max_labels of them keep the top max_labels by supported k-mers, ties by column id, then
    seq_id, and state labels_truncated). The source record's own label supports every in-graph k-mer of a query taken
    from it (the graph_runs summed; under a stop, the resolved prefix's), so a truncated list may leave it out only
    in a tie at that support: min_kept_kmers and max_dropped_kmers equal to it and dropped_full_length at least 1.
    Such an answer passes, the ranking being all an answer can show (the tie's column ids are not in it); a list that
    keeps a label with fewer k-mers, ranks a dropped label above a kept one or drops none at full support fails, as
    does a complete list without the label. Only more than max_labels labels are truncated, to exactly max_labels
    (the request's, else the default): a list cut to another size fails, tie or not. A discovery stopped before it
    read the row of an in-graph k-mer (SPEC-labeled-traversal-core §4.5: the resolve of a prefix without one) has no
    profile to rank: inconclusive."""
    body = _body(ans)
    names = [l.get("label") for l in _get(body, "labels")]
    trunc = body.get("labels_truncated")
    if a["label"] in names:
        return PASS, "%s discovered (%d labels%s)" % (a["label"], len(names), ", truncated" if trunc else "")
    if body.get("stop") and not trunc and not names and _in_graph_kmers(body) == 0:
        return INCONCLUSIVE, "%s not listed: the discovery stopped before it read an in-graph k-mer (resolved_kmers " \
                             "%r), so it has no profile to rank" % (
                                 a["label"], (body["stop"] or {}).get("resolved_kmers"))
    if not trunc:
        return FAIL, "%s not discovered in a complete list of %d labels" % (a["label"], len(names))
    if not isinstance(trunc, dict):
        return FAIL, "labels_truncated %r, not the SPEC's object" % (trunc,)
    vals = {}
    for f in ("kept", "total", "min_kept_kmers", "max_dropped_kmers", "dropped_full_length"):
        if f not in trunc:
            raise Missing("labels_truncated." + f)
        vals[f] = _number(trunc[f])
    full = _in_graph_kmers(body)
    stated = "kept %s of %s, min_kept_kmers %s, max_dropped_kmers %s, dropped_full_length %s" % (
        trunc["kept"], trunc["total"], trunc["min_kept_kmers"], trunc["max_dropped_kmers"],
        trunc["dropped_full_length"])
    if None in vals.values():
        return FAIL, "%s not listed, labels_truncated not numeric (%s)" % (a["label"], stated)
    want = _discover_max_labels(sc.get("request") if isinstance(sc.get("request"), dict) else {})
    if want is not None and (vals["kept"] != want or vals["total"] <= want):
        return FAIL, "%s not listed, and the list is truncated to %s of %s labels for max_labels %d: only more than " \
                     "max_labels labels are truncated, to exactly max_labels (SPEC-labeled-traversal-core §4.2)" % (
                         a["label"], trunc["kept"], trunc["total"], want)
    if vals["max_dropped_kmers"] > vals["min_kept_kmers"]:
        return FAIL, "%s not listed, and the truncation is not by support: a dropped label has more supported " \
                     "k-mers than a kept one (%s; SPEC-labeled-traversal-core §4.2)" % (a["label"], stated)
    if vals["min_kept_kmers"] < full:
        return FAIL, "%s not listed, but it supports all %d in-graph k-mers and a kept label only %s (%s): the " \
                     "top max_labels by supported k-mers hold it (SPEC-labeled-traversal-core §4.2)" % (
                         a["label"], full, trunc["min_kept_kmers"], stated)
    if vals["dropped_full_length"] < 1 or vals["max_dropped_kmers"] < full:
        return FAIL, "%s not listed, but no dropped label supports all %d in-graph k-mers as it does (%s)" % (
            a["label"], full, stated)
    return PASS, "%s not listed: it ties the kept labels at all %d in-graph k-mers (%s), a tie decided by column " \
                 "id, then seq_id (SPEC-labeled-traversal-core §4.2), which the answer does not state" % (
                     a["label"], full, stated)


def _discovery_problems(req: dict, body: dict) -> List[str]:
    """A discovery's list and labels_truncated as SPEC-labeled-traversal-core §4.2 ranks them: the labels by
    supported k-mers, most first (each label's kmers_supported is that count unless support is trace, whose count
    is of trace runs), at most max_labels (the request's, else the default); a truncated list states kept (its length,
    exactly max_labels) of total (more), the least support kept at least the most dropped, and full-length drops
    only where a dropped label has every in-graph k-mer, none of which is more than the query has."""
    probs: List[str] = []
    labels = body.get("labels") or []
    trunc = body.get("labels_truncated")
    want = _discover_max_labels(req)
    full = _in_graph_kmers(body)
    ranked = req.get("support") != "trace"
    ks = [_number(l.get("kmers_supported")) for l in labels]
    if want is not None and len(labels) > want:
        probs.append("%d labels discovered, max_labels %d" % (len(labels), want))
    if ranked and None not in ks and any(x < y for x, y in zip(ks, ks[1:])):
        probs.append("discovered labels not ranked by supported k-mers: %s" % ks[:8])
    if any(x is not None and x > full for x in ks):
        probs.append("a discovered label supports more k-mers than the %d in the graph" % full)
    if trunc is None:
        return probs
    if not isinstance(trunc, dict):
        return probs + ["labels_truncated %r, not the SPEC's object" % (trunc,)]
    v = {f: _number(trunc.get(f)) for f in ("kept", "total", "min_kept_kmers", "max_dropped_kmers",
                                            "dropped_full_length")}
    if None in v.values():
        return probs + ["labels_truncated %s, a field missing or not numeric" % json.dumps(trunc, sort_keys=True)]
    if v["kept"] != len(labels) or v["total"] <= v["kept"]:
        probs.append("labels_truncated kept %s of total %s for %d labels listed" % (
            trunc["kept"], trunc["total"], len(labels)))
    if want is not None and (v["kept"] != want or v["total"] <= want):
        probs.append("labels_truncated kept %s of total %s with max_labels %d: a discovery truncates only more than "
                     "max_labels labels, to exactly max_labels (SPEC-labeled-traversal-core §4.2)" % (
                         trunc["kept"], trunc["total"], want))
    if v["max_dropped_kmers"] > v["min_kept_kmers"]:
        probs.append("labels_truncated max_dropped_kmers %s above min_kept_kmers %s: not truncated by support" % (
            trunc["max_dropped_kmers"], trunc["min_kept_kmers"]))
    if v["min_kept_kmers"] > full:
        probs.append("labels_truncated min_kept_kmers %s above the %d in-graph k-mers" % (
            trunc["min_kept_kmers"], full))
    if not 0 <= v["dropped_full_length"] <= v["total"] - v["kept"] or (
            v["dropped_full_length"] > 0 and v["max_dropped_kmers"] != full):
        probs.append("labels_truncated dropped_full_length %s with max_dropped_kmers %s, %d in-graph k-mers, %s "
                     "dropped" % (trunc["dropped_full_length"], trunc["max_dropped_kmers"], full,
                                  v["total"] - v["kept"]))
    if ranked and labels and None not in ks and ks[-1] != v["min_kept_kmers"]:
        probs.append("labels_truncated min_kept_kmers %s, the last kept label supports %s" % (
            trunc["min_kept_kmers"], ks[-1]))
    return probs


def _resolve_invariants(sc, ans, a, arm, client=None):
    if ans.get("http") != 200:
        return PASS, "no answer to check (http %s)" % ans.get("http")
    body = _body(ans)
    req = sc["request"]
    probs = []
    stop = body.get("stop")
    seq = req.get("sequence", "") if isinstance(req, dict) else ""
    qk = max(0, len(seq) - K + 1)
    num = body.get("num_kmers")
    # the deadline's terms: a clamp only of what the route clamps, no budget lowered unstated (the stop rule reads
    # limits.time_budget_ms), and the host's reserve (L1 of the verification of SD-021)
    probs += _resolve_clamp_problems(sc, body, arm)
    reserve = _resolve_reserve_problem(body)
    if reserve:
        probs.append(reserve)
    if stop:
        bad = _resolve_stop_problem(sc, body, client)
        if bad and not (reserve and bad.endswith(reserve)):
            probs.append(bad)
        x = stop.get("resolved_kmers")
        if num != x:
            probs.append("num_kmers %s != stop.resolved_kmers %s" % (num, x))
        if stop.get("remainder_from_bp") != x:
            probs.append("remainder_from_bp %s != resolved_kmers %s" % (stop.get("remainder_from_bp"), x))
        if stop.get("resolved_bp") != (x + K - 1 if x else 0):
            probs.append("resolved_bp %s for resolved_kmers %s" % (stop.get("resolved_bp"), x))
        if stop.get("query_kmers") != qk or stop.get("query_bp") != len(seq):
            probs.append("query_kmers/bp %s/%s" % (stop.get("query_kmers"), stop.get("query_bp")))
    elif num != qk:
        probs.append("num_kmers %s, expected %d" % (num, qk))
    gr = body.get("graph_runs") or []
    prev = -1
    for r in gr:
        if not (r[0] < r[1] and r[0] > prev and r[1] <= (num or 0)):
            probs.append("graph_runs not sorted/disjoint within [0, num_kmers): %s" % gr[:6])
            break
        prev = r[1]
    for lab in body.get("labels") or []:
        runs = lab.get("runs") or []
        if isinstance(runs, list) and runs and isinstance(runs[0], list):
            total = sum(r[1] - r[0] for r in runs)
            if total != lab.get("kmers_supported"):
                probs.append("label %s: runs sum %d != kmers_supported %s" % (lab.get("label"), total, lab.get("kmers_supported")))
            for r in runs:
                if not any(g[0] <= r[0] and r[1] <= g[1] for g in gr):
                    probs.append("label %s: run %s outside the graph runs" % (lab.get("label"), r))
                    break
    if isinstance(req, dict) and isinstance(req.get("discover"), dict):
        probs += _discovery_problems(req, body)
    if isinstance(req, dict) and req.get("labels") and not probs:
        names = [l.get("label") for l in body.get("labels") or []]
        missing = [l for l in req["labels"] if l not in names]
        if missing:
            probs.append("requested labels missing from the answer: %s" % missing[:5])
    return (FAIL, "; ".join(probs[:6])) if probs else (PASS, "resolve answer consistent")


# ---------------------------------------------------------------------------------------------- traverse checks
_END_CODES = {"D": "dead_end", "L": "label_lost", "B": "loss_budget", "R": "branch", "U": "edge_reuse",
              "V": "edge_reuse_rc", "J": "rejoined_seed", "T": "trace_break", "X": "max_extension_bp",
              "S": "max_steps", "P": "max_live_paths", "N": "max_paths", "O": "max_output_bp", "M": "time_budget",
              "W": "beam_pruned", "Y": "resource_limit"}
_END_QUALIFIERS = {"Rm": "minority", "Rb": "below_min_labels", "Rs": "split_limit", "Dh": "hairpin",
                   "Ls": "superseded", "Lx": "switch_sources", "Lw": "label_lost"}


_MGT_INT = re.compile(r"(0|[1-9][0-9]*)\Z")
_MGT_RANGES_ITEM = re.compile(r"(0|[1-9][0-9]*)(?:-(0|[1-9][0-9]*))?\Z")
_MGT_U64_DIGITS, _MGT_U64 = 20, (1 << 64) - 1
_MGT_UNPCT = {"%25": "%", "%0A": "\n", "%0D": "\r"}


def _mgt_ranges(tok: str, bound: Optional[int]) -> set:
    """An MGT RANGES token (SPEC §7.5.1: one valid spelling, every other one a ValueError): `.` empty, else ascending
    collapsed items `a` or `a-b` (b > a, no leading zero), each id 64-bit and below |bound| (the number of L
    records), checked before a run is expanded (with no bound, a run of more than 10**7 ids is refused). The library's
    _codec.decode_ranges, held to the golden vectors by tests/test_checks.py MgtCodecVectors."""
    if tok == ".":
        return set()
    out: set = set()
    last = None
    for item in tok.split(","):
        m = _MGT_RANGES_ITEM.match(item)
        if not m or max(len(m.group(1)), len(m.group(2) or "")) > _MGT_U64_DIGITS:
            raise ValueError("RANGES item %r: an id is a decimal integer without sign or leading zero, a run a-b, the "
                             "items comma-separated (one spelling per token)" % item)
        lo = int(m.group(1))
        hi = int(m.group(2)) if m.group(2) is not None else lo
        if m.group(2) is not None and hi <= lo:
            raise ValueError("RANGES run %r: a run a-b has b above a" % item)
        if last is not None and lo <= last + 1:
            raise ValueError("RANGES: ids strictly ascending, consecutive ids collapsed into a run a-b")
        if hi > _MGT_U64:
            raise ValueError("RANGES item %r: an id is a 64-bit integer" % item)
        if bound is not None and hi >= bound:
            raise ValueError("RANGES item %r: every id below the number of L records (%d)" % (item, bound))
        if bound is None and hi - lo > 10 ** 7:
            raise ValueError("RANGES run %r: more than 10**7 ids with no id space to check them against" % item)
        out.update(range(lo, hi + 1))
        last = hi
    return out


def _mgt_setexpr(tok: str, base: Optional[set], bound: Optional[int]) -> set:
    """An MGT SETEXPR (SPEC §7.5.1): RANGES, or a delta `!`, `!R`, `!+A`, `!R+A` against |base| (an empty part is
    omitted, never `.`); a delta whose R is not inside the base or whose A meets it is a ValueError: the reader
    derived another base than the writer (the library's _codec.decode_setexpr)."""
    if not tok.startswith("!"):
        return _mgt_ranges(tok, bound)
    if base is None:
        raise ValueError("SETEXPR: a delta needs a base (a root's entry set is explicit only)")
    removed, plus, added = tok[1:].partition("+")
    if removed == "." or (plus and added in ("", ".")):
        raise ValueError("SETEXPR: an empty part of a delta is omitted, never written '.' or left after '+'")
    r = _mgt_ranges(removed, bound) if removed else set()
    a = _mgt_ranges(added, bound) if plus else set()
    if not r <= base or a & base:
        raise ValueError("SETEXPR: a delta's R lies in its base and its A outside it (else the reader derived "
                         "another base)")
    return (set(base) - r) | a


def _mgt_unescape(text: str) -> str:
    """A free-text last field (SPEC §7.5.1): %25, %0A and %0D are its only escapes; any other `%` sequence or a raw
    LF/CR is a ValueError (the library's _codec.pct_unescape)."""
    if "\n" in text or "\r" in text:
        raise ValueError("free text: no raw LF or CR (written %0A, %0D)")
    if "%" not in text:
        return text
    out, i = [], 0
    while i < len(text):
        if text[i] == "%":
            esc = text[i:i + 3]
            if esc not in _MGT_UNPCT:
                raise ValueError("free text: the escape %r; %%25, %%0A and %%0D are the only escapes" % esc)
            out.append(_MGT_UNPCT[esc])
            i += 3
        else:
            out.append(text[i])
            i += 1
    return "".join(out)


_MGT_FRONT = re.compile(r"(0|[1-9][0-9]*) (.*)\Z", re.S)


def _mgt_front(prev: bytes, tail: str) -> bytes:
    """An L record's name from its tail `<prefix_len> <suffix>` and the previous name's UTF-8 |prev| (SPEC §7.5.2):
    the previous name's first prefix_len bytes and the unescaped suffix; a tail without its separating space, a
    prefix that is no integer, goes beyond the previous name or ends inside a character is a ValueError (the
    library's _codec.front_decode_parts)."""
    m = _MGT_FRONT.match(tail)
    if not m or len(m.group(1)) > _MGT_U64_DIGITS:
        raise ValueError("front coding: '<prefix_len> <suffix>', prefix_len an integer without sign or leading zero")
    p = int(m.group(1))
    if p > len(prev):
        raise ValueError("front coding: prefix_len %d beyond the previous name's %d bytes" % (p, len(prev)))
    name = prev[:p] + _mgt_unescape(m.group(2)).encode("utf-8")
    try:
        name.decode("utf-8")
    except UnicodeDecodeError:
        raise ValueError("front coding: prefix_len %d ends inside a character of the previous name" % p) from None
    return name


def graphlet_label_dict(mgt: str) -> List[dict]:
    """The label dictionary of an MGT v1 body, from its L records (SPEC §7.5.2: name = the previous name's first
    prefix_len UTF-8 bytes + the suffix; label id = ordinal): [{name, kind header|column, column, seq_id}]. A record
    the strict reader rejects is a ValueError naming the record, its token and the rule (SD-014)."""
    out: List[dict] = []
    prev = b""
    for line in mgt.split("\n"):
        if not line.startswith("L "):
            continue
        f = line.split(" ", 4)
        where = "the L record of label id %d" % len(out)
        if len(f) < 5:
            raise ValueError("%s: %d fields, not L <c|h> <column> <seq_id|*> <prefix_len> <suffix>" % (where, len(f)))
        try:
            name = _mgt_front(prev, f[4])
        except ValueError as e:
            raise ValueError("%s, tail %r: %s" % (where, f[4], e)) from None
        if f[3] != "*" and not _MGT_INT.match(f[3]):
            raise ValueError("%s, seq_id %r: an integer without sign or leading zero, or *" % (where, f[3]))
        prev = name
        out.append({"name": name.decode("utf-8"), "kind": "header" if f[1] == "h" else "column",
                    "column": int(f[2]) if f[2].isdigit() else f[2], "seq_id": None if f[3] == "*" else int(f[3])})
    return out


def _mgt_token(fn: Callable[..., set], where: str, field: str, tok: str, *args: Any) -> set:
    """fn(tok, *args) of a label token of record |where|; a ValueError restated with the record, the field and the
    token (SD-014: the check that reads the sets names all three)."""
    try:
        return fn(tok, *args)
    except ValueError as e:
        raise ValueError("%s, %s %r: %s" % (where, field, tok, e)) from None


def _graphlet_label_sets(a: dict, mode: str, num_seed_labels: Optional[int], bound: Optional[int]) -> None:
    """Each segment's entry set as `labels` (and labels_total, labels_truncated) of a `detail: full` segment (SPEC
    §7.1: the names live at the segment's first step), derived from the G records as §7.5.3 states it (the
    reference: the library's derive.py): G.entry against its base (parent[0]'s end set, a merge's the union of its
    parents' end sets, the root none), `*` the seed labels at the root and a merge's partition union (constrain);
    G.end against the entry, `*` the entry with, by position and run ends before switch-ins, every run anchored on
    the segment ending before its last node removed and every switch target added (constrain), the last P set or the
    entry (annotate). Raises ValueError on a body the rules do not read, naming the record, its token and the rule."""
    segs = a["segments"]
    if mode == "c" and num_seed_labels is not None and bound is not None and num_seed_labels > bound:
        raise ValueError("the S record, num_seed_labels %d: at most the number of L records (%d)" % (
            num_seed_labels, bound))
    runs_on: Dict[int, List[dict]] = {}
    for r in a["_runs"]:
        runs_on.setdefault(r["segment"], []).append(r)
    end: Dict[int, set] = {}
    for sg in segs:
        g = sg.pop("_g")
        parents = sg["parents"]
        where = "the G record of %s segment %d" % (a.get("_side"), sg["id"])
        if any(p >= sg["id"] or p not in end for p in parents):
            raise ValueError("%s, parents %s: every parent an earlier segment" % (where, parents))
        base = None if not parents else set().union(*(end[p] for p in parents))
        if g["entry"] == "*":
            if not parents and mode == "c" and num_seed_labels is not None:
                entry = set(range(num_seed_labels))
            elif len(parents) > 1 and mode == "c" and g["partition"] != "*":
                entry = set().union(*(_mgt_token(_mgt_ranges, where, "partition", x, bound)
                                      for x in g["partition"].split("|")))
            else:
                raise ValueError("%s, entry '*': '*' derives the root's entry from the S record's num_seed_labels and "
                                 "a merge's from its partition, in constrain mode only (SPEC §7.5.3)" % where)
        else:
            entry = _mgt_token(_mgt_setexpr, where, "entry", g["entry"], base, bound)
        if g["entry_total"] != "*" and not _MGT_INT.match(g["entry_total"]):
            raise ValueError("%s, entry_total %r: an integer without sign or leading zero, or *" % (
                where, g["entry_total"]))
        total = len(entry) if g["entry_total"] == "*" else int(g["entry_total"])
        if g["end"] != "*":
            e = _mgt_token(_mgt_setexpr, where, "end", g["end"], entry, bound)
        elif mode == "c":
            stop = sg["from_bp"] + sg["length_bp"]
            ops = [(r["to_bp"], 0, r["label"]) for r in runs_on.get(sg["id"], []) if r["to_bp"] < stop]
            ops += [(at, 1, to) for at, to in g["switches"]]
            e = set(entry)
            for _, kind, label in sorted(ops, key=lambda o: (o[0], o[1])):
                if kind == 0:
                    e.discard(label)
                else:
                    e.add(label)
        else:
            e = set(entry)
            for j, ptok in enumerate(g["presence"]):
                e = _mgt_token(_mgt_setexpr, "the P record %d of %s segment %d" % (j, a.get("_side"), sg["id"]),
                               "SETEXPR", ptok, e, bound)
        end[sg["id"]] = e
        sg["labels"] = sorted(entry)
        sg["labels_total"] = total
        sg["labels_truncated"] = total > len(entry)


def graphlet_arms(mgt: str) -> Dict[str, dict]:
    """The arms of an MGT v1 body (SPEC-labeled-traversal-core.md §7.5) in the shape of a `detail: full` arm, as
    far as the checks read it: segments (id, parents, from_bp, length_bp, sequence in natural orientation, events:
    the blocked successors of the E records and the label_end events of the R records with their structural
    successors; labels, labels_total and labels_truncated: the entry sets of the G records, _graphlet_label_sets,
    left out when the body's sets cannot be read, the arm's _labels_rejected then saying why), splits (non-empty when a segment ends in a split), branch_events
    (the V records), and paths (one per leaf, the chain through parents[0], length_bp = from_bp + length_bp,
    end_reasons and end_labels (label, branches) from the runs anchored at the leaf with to_bp = the path length,
    path_reason from T), and the A record's counters and cap_trigger (reason, at_bp)."""
    arms: Dict[str, dict] = {}
    cur = None
    seg = None
    mode, num_seed_labels, n_labels = "c", None, 0
    for line in mgt.split("\n"):
        if not line:
            continue
        f = line.split(" ")
        rec = f[0]
        if rec == "H" and len(f) > 6:
            mode = f[6]
        elif rec == "S" and len(f) > 4 and f[4].isdigit():
            num_seed_labels = int(f[4])
        elif rec == "L":
            n_labels += 1
        if rec == "A":
            side = "left" if f[1] == "l" else "right"
            cur = {"status": {"c": "complete", "t": "truncated", "p": "pruned"}.get(f[2], f[2]),
                   "segments": [], "splits": [], "branch_events": [], "paths": [], "_runs": [], "_leaves": {}}
            # the arm's counters (name=value,...) and cap trigger (reason,at_bp,segment,live_paths,live_labels,
            # exact,demand, or *) of SPEC §7.5's A record; a result's JSON arm states them too and is preferred
            if len(f) > 11:
                cnt = {}
                for kv in f[10].split(","):
                    name, _, val = kv.partition("=")
                    if name and val.lstrip("-").isdigit():
                        cnt[name] = int(val)
                if cnt:
                    cur["counters"] = cnt
                if f[11] != "*":
                    c = f[11].split(",")
                    cur["cap_trigger"] = {"reason": _END_CODES.get(c[0], c[0]),
                                          "at_bp": int(c[1]) if len(c) > 1 and c[1].isdigit() else None}
            arms[side] = cur
            cur["_side"] = side
            seg = None
        elif cur is None:
            continue
        elif rec == "G":
            if len(f) != 11:
                raise Missing("G record of %d fields" % len(f))
            parents = [] if f[1] == "*" else [int(x) for x in f[1].split(",")]
            bases = "" if f[10] in (".", "*") else f[10]
            natural = bases if cur["_side"] == "right" else bases[::-1]
            # every segment lists its events (the E-b and R records below add them): an MGT states all of them
            # (SPEC §7.5), so a segment without any has none, as a JSON tree's empty list says
            seg = {"id": len(cur["segments"]), "parents": parents, "from_bp": int(f[2]), "length_bp": int(f[3]),
                   "sequence": natural, "events": [],
                   "_g": {"entry": f[4], "entry_total": f[5], "end": f[6], "partition": f[7], "switches": [],
                          "presence": []}}
            if f[10] == "*":
                del seg["sequence"]
            cur["segments"].append(seg)
            if f[8] in ("0", "1"):
                cur["splits"].append({"segment": seg["id"], "ambiguous": f[8] == "1"})
        elif rec == "V":
            cur["branch_events"].append({"at_bp": int(f[1]), "segment": int(f[2])})
        elif rec == "E" and seg is not None and len(f) >= 5 and f[2] == "s":
            seg["_g"]["switches"].append((int(f[1]), int(f[4])))   # a switch: E <at_bp> s <from> <to> <cost>
        elif rec == "P" and seg is not None and len(f) >= 5:
            seg["_g"]["presence"].append(f[4])   # an annotate run: P <from_bp> <to_bp> <total> <SETEXPR>
        elif rec == "E" and seg is not None and len(f) >= 7 and f[2] == "h":
            # a hairpin step: E <at_bp> h <char> <total> <RANGES> <f|s>, followed or skipped (SPEC §7.5.2), the
            # JSON event's at_bp, char and followed
            seg.setdefault("events", []).append({"at_bp": int(f[1]), "type": "hairpin", "char": f[3],
                                                 "followed": f[6] == "f"})
        elif rec == "E" and seg is not None and len(f) >= 5 and f[2] == "b":
            # a blocked successor: E <at_bp> b <char> <reason code> ... (SPEC §7.5)
            seg.setdefault("events", []).append({"at_bp": int(f[1]), "type": "blocked", "char": f[3],
                                                 "reason": _END_CODES.get(f[4][0], f[4])})
        elif rec == "T" and seg is not None:
            cur["_leaves"][seg["id"]] = None if f[1] == "*" else _END_CODES.get(f[1][0], f[1])
        elif rec == "R":
            run = {"segment": int(f[1]), "label": int(f[2]), "from_bp": int(f[3]), "to_bp": int(f[4]), "end": f[5]}
            # R ... <structural_successors|*> <branches> <loss>: the lineage's values when the run ended (§7.5.2)
            if len(f) > 10 and f[9].isdigit():
                run["structural_successors"] = int(f[9])
            if len(f) > 10 and f[10].isdigit():
                run["branches"] = int(f[10])
            cur["_runs"].append(run)
    for side, a in arms.items():
        segs = a["segments"]
        try:
            _graphlet_label_sets(a, mode, num_seed_labels, n_labels)
        except (ValueError, KeyError, IndexError, TypeError) as e:
            # a body the strict reader rejects (SPEC §7.5.1, §7.5.3): no label sets, and the checks that read them
            # fail naming why (_mgt_rejected; SD-014: never the sequence alone)
            a["_labels_rejected"] = str(e) if isinstance(e, ValueError) else "%s: %s" % (type(e).__name__, e)
            for sg in segs:
                for key in ("_g", "labels", "labels_total", "labels_truncated"):
                    sg.pop(key, None)
        # the label_end events (SPEC §7.5.6: derived from the R records with a code end, on their segment at to_bp,
        # reason = the qualifier's text or the code's enum; `Lw` ends silently and `m` is closed by a merge)
        for r in a["_runs"]:
            end = r["end"]
            if end in ("m", "Lw") or not 0 <= r["segment"] < len(segs):
                continue
            ev = {"at_bp": r["to_bp"], "type": "label_end", "label": r["label"],
                  "reason": _END_QUALIFIERS.get(end) or _END_CODES.get(end[0], end)}
            if "structural_successors" in r:
                ev["structural_successors"] = r["structural_successors"]
            segs[r["segment"]].setdefault("events", []).append(ev)
        for leaf, reason in sorted(a["_leaves"].items()):
            chain = [leaf]
            while segs[chain[-1]]["parents"]:
                chain.append(segs[chain[-1]]["parents"][0])
            chain.reverse()
            length = segs[leaf]["from_bp"] + segs[leaf]["length_bp"]
            ends: Dict[str, int] = {}
            end_labels = []
            for r in a["_runs"]:
                if r["segment"] == leaf and r["to_bp"] == length and r["end"] != "m":
                    # the label's end reason as the full result counts it: the code's enum (Dh -> dead_end)
                    name = _END_CODES.get(r["end"][0], r["end"])
                    ends[name] = ends.get(name, 0) + 1
                    end_labels.append({"label": r["label"], "branches": r.get("branches")})
            a["paths"].append({"id": len(a["paths"]), "segments": chain, "length_bp": length, "end_reasons": ends,
                               "end_labels": sorted(end_labels, key=lambda x: x["label"]), "path_reason": reason})
        for k in ("_runs", "_leaves", "_side"):
            a.pop(k, None)
    return arms


def _result(ans: dict, i: int) -> dict:
    """The i-th seed result; a graphlet result (`detail: graphlet`, or the MCP's graphlet export) is read through
    its MGT body into the arm shape of `detail: full`."""
    res = _get(_body(ans), "results")
    if i >= len(res):
        raise Missing("results[%d]" % i)
    r = res[i]
    if isinstance(r.get("graphlet"), str) and not any("segments" in a for a in (r.get("arms") or {}).values()):
        r = dict(r)
        parsed = graphlet_arms(r["graphlet"])
        rejected = []   # why the strict reader rejects the body's label dictionary or an arm's label sets
        try:
            ld = graphlet_label_dict(r["graphlet"])
            if "label_dict" not in r:
                r["_mgt_label_dict"] = ld
        except ValueError as e:
            rejected.append(str(e))
        for side in ("left", "right"):
            why = (parsed.get(side) or {}).get("_labels_rejected")
            if why and why not in rejected:   # a body-level rule (the S record's) rejects both arms alike
                rejected.append(why)
        if rejected:
            r["_mgt_rejected"] = "; ".join(rejected)
        arms = {}
        for side, a in parsed.items():
            merged = dict((r.get("arms") or {}).get(side) or {})
            for key, value in a.items():
                if key in ("counters", "cap_trigger") and key in merged:
                    continue   # the JSON arm's own statement (exact types) over the MGT's copy
                merged[key] = value
            arms[side] = merged
        r["arms"] = arms
    return r


def _arm_flanks(arm_obj: dict, side: str) -> List[Tuple[dict, str]]:
    """(path, flank) for each path of an arm, the flank in natural orientation (SPEC §7.1: a right flank is the
    concatenation root -> leaf, a left flank leaf -> root)."""
    segs = {s["id"]: s for s in arm_obj.get("segments", [])}
    out = []
    for p in arm_obj.get("paths", []):
        parts = []
        for sid in p.get("segments", []):
            s = segs.get(sid)
            if s is None or "sequence" not in s:
                raise Missing("segment sequences")
            parts.append(s["sequence"])
        flank = "".join(parts) if side == "right" else "".join(reversed(parts))
        out.append((p, flank))
    return out


def _seg_events(seg: Any) -> List[dict]:
    """The events of a segment that a check can read: its `events` list's objects; none where the segment states no
    list (absent, null or another type: _structural_end_problems requires the list where an end needs it)."""
    evs = seg.get("events") if isinstance(seg, dict) else None
    return [x for x in evs if isinstance(x, dict)] if isinstance(evs, list) else []


def _leaf_ends(armo: dict, p: dict, path_reason: bool = True) -> set:
    """The end reasons of path |p|'s leaf: its path_reason (unless |path_reason| is False), its end_reasons and the
    reasons of the label_end events on its leaf segment at the path's end (SPEC-labeled-traversal-core §7.1, §7.5.6:
    a run that ends at the leaf states its label_end event there at to_bp = the path length; a quorum or split-limit
    stop is the end reason `branch` with the qualifier text, e.g. below_min_labels, in that event: §13, I3 note 2).
    Every route to a leaf ends with the leaf's reasons, whichever parent it took at a merge."""
    ends = set(p.get("end_reasons") or {})
    if path_reason:
        ends.add(p.get("path_reason"))
    segs = p.get("segments") or []
    if segs:
        leaf = next((x for x in armo.get("segments") or [] if x.get("id") == segs[-1]), None)
        for ev in _seg_events(leaf):
            if ev.get("type") == "label_end" and ev.get("reason") and ev.get("at_bp") == p.get("length_bp"):
                ends.add(ev["reason"])
    ends.discard(None)
    return ends


# A merged walk (branching.on_reconverge "merge", the walk preset's default) is displayed through the first parent of
# each merge (SPEC-labeled-traversal-core §7.1: the parent carried by the most labels, or the arrival order under a
# memory budget in tree / full detail), so after a split that several records (or repeat copies of one) carry, the
# displayed path may be another record's allele while the record's own route passes the merge through another parent.
# The walk checks therefore search the arm's segment DAG through segments[].parents, taking ANY parent at a merge
# (labels_via_parent is not followed: a label's kept entry may come through the other allele's parent). Where the
# segments state their labels (segments[].labels of `detail: full`, the G records' entry sets of a graphlet) a route
# passes only segments that carry the record's own label (its header or its column, by the label dictionary; a
# segment whose list the cap cut, labels_truncated, passes): the record's sequence must not be pieced together from
# other labels' segments (A-04). The search is pruned at the first base that departs from the record and memoised on
# (segment, the seed copies still followed): every route to a segment has the same length, and the copies it still
# follows fix the spelled flank, so each such state is searched once (A-05). _MAX_ROUTE_STATES bounds the states as a
# safeguard; a check that could not finish the search is inconclusive, never failed.
_MAX_ROUTE_STATES = 100000


def _mgt_rejected(r: dict) -> Optional[str]:
    """Why the strict MGT reader rejects the body of graphlet result |r| (_result: its L records, or an arm's label
    tokens), as the detail of the check that reads its labels, or None. The owner's decision on SD-014: such a body
    fails walk_in_record and walk_end on every arm, naming the record, the token and the rule; never the sequence
    alone (the fallback before)."""
    why = r.get("_mgt_rejected")
    if not why:
        return None
    return "the graphlet body breaks SPEC-labeled-traversal-core §7.5: %s (its labels are not read, and the walk is " \
           "not judged on the sequence alone)" % why


def _record_label_ids(r: dict, t: dict) -> Optional[set]:
    """The label ids of the truth's record in a seed result's label dictionary (its header, or its taxid's column),
    the graphlet's L records where the result has no label_dict; None without any (the label rule then does not
    apply, see _route_labels)."""
    ld = r.get("label_dict")
    if not isinstance(ld, list):
        ld = r.get("_mgt_label_dict")
    if not isinstance(ld, list):
        return None
    names = set(str(t[k]) for k in ("record", "taxid") if t.get(k) is not None)
    return set(i for i, x in enumerate(ld) if isinstance(x, dict) and str(x.get("name")) in names)


def _route_labels(r: dict, t: dict, armo: dict) -> Tuple[Optional[set], str]:
    """(the record's label ids the route search requires, or None for the sequence-only rule, and why it is
    sequence-only: '' when the label rule applies)."""
    ids = _record_label_ids(r, t)
    if not ids:
        return None, "sequence only: the record's label is not in the seed's label dictionary"
    if not any(isinstance(x.get("labels"), list) for x in armo.get("segments") or []):
        return None, "sequence only: the segments state no labels"
    return ids, ""


def _record_routes(armo: dict, side: str, frames: Sequence[Tuple[Callable[[int], str], int, int]],
                   cap: Optional[int] = None, labels: Optional[set] = None,
                   event: Optional[Tuple[int, str, str]] = None) -> Tuple[List[Tuple[dict, str]], int, bool]:
    """The root-to-leaf routes of arm |armo| through segments[].parents, any parent at a merge, whose spelled flank
    (natural orientation: a right flank root -> leaf, a left flank leaf -> root) is, for its whole length, the
    continuation of one of |frames| (_seed_frame: one per copy of the seed; (cont, room, end)). With |labels| (label
    ids), a segment that lists its labels without the cap having cut them must carry one of them. With |event|
    (at_bp, char, reason), a route counts only once one of its segments states that blocked event. Returns (routes as
    (the leaf's paths[] entry, flank), the number of routes that followed a frame to the end of its truth flank and
    ran past it, whether the search stopped at |cap| states, default _MAX_ROUTE_STATES)."""
    cap = _MAX_ROUTE_STATES if cap is None else cap
    segs = dict((x.get("id"), x) for x in armo.get("segments") or [])
    children: Dict[Any, List[Any]] = {}
    for x in armo.get("segments") or []:
        for par in x.get("parents") or []:
            children.setdefault(par, []).append(x.get("id"))
    leaves = dict((p["segments"][-1], p) for p in armo.get("paths") or [] if p.get("segments"))
    stack = [(x.get("id"), "", tuple(range(len(frames))), event is None)
             for x in reversed(armo.get("segments") or []) if not x.get("parents")]
    routes: List[Tuple[dict, str]] = []
    beyond = visits = 0
    seen = set()
    while stack:
        sid, flank, alive, stated = stack.pop()
        if (sid, alive, stated) in seen:
            continue   # the same segment, copies and statement: the same flank and subtree (A-05)
        seen.add((sid, alive, stated))
        visits += 1
        if visits > cap:
            return routes, beyond, True
        sg = segs.get(sid)
        if sg is None or "sequence" not in sg:
            raise Missing("segment sequences")
        if labels is not None and isinstance(sg.get("labels"), list) and not sg.get("labels_truncated") and \
                not labels.intersection(sg["labels"]):
            continue
        if not stated:
            stated = any(ev.get("type") == "blocked" and ev.get("at_bp") == event[0] and ev.get("char") == event[1]
                         and ev.get("reason") == event[2] for ev in _seg_events(sg))
        flank = flank + sg["sequence"] if side == "right" else sg["sequence"] + flank
        n = len(flank)
        keep, past = [], False
        for fi in alive:
            cont, room, _ = frames[fi]
            if room < n:
                # the part inside the truth flank follows it: the route runs past what the truth can judge
                inside = flank[:room] if side == "right" else flank[n - room:]
                past = past or inside == cont(room)
            elif flank == cont(n):
                keep.append(fi)
        if not keep:
            beyond += 1 if past and stated else 0
            continue
        kids = children.get(sid) or []
        if not kids and stated:
            routes.append((leaves.get(sid) or {"id": None, "segments": [sid], "length_bp": n}, flank))
        for c in reversed(kids):
            stack.append((c, flank, tuple(keep), stated))
    return routes, beyond, False


def _arms_of(sc: dict, i: int) -> List[str]:
    d = (sc["request"].get("strategy") or {}).get("direction", "both")
    return ["left", "right"] if d == "both" else [d]


# path ends that may legitimately stop a walk along its record short of the radius and the record's end: a stated
# branch (a split the strategy does not follow further) or the walk meeting its own edges again (a repeat), each
# only where it can have happened (_structural_end_problems)
_STATED_STRUCTURAL_ENDS = {"branch", "edge_reuse", "edge_reuse_rc", "rejoined_seed", "minority", "below_min_labels",
                           "split_limit"}
# a stated resource limit: the walk stopped before it could show the property (inconclusive, never a pass)
_RESOURCE_ENDS = {"max_steps", "max_live_paths", "max_paths", "max_output_bp", "time_budget", "resource_limit",
                  "beam_pruned"}
# the quorum and split-limit ends of SPEC-labeled-traversal-core §6.4: a path counts them as `branch`
# (end_reasons; §13 I3 note 2), its label_end event names the stop (`reason`, or `text` beside reason `branch`; a
# graphlet's R qualifier Rm / Rb / Rs)
_QUORUM_ENDS = {"minority", "below_min_labels", "split_limit"}
# the ends of a step not taken (§6.6), each stated by a blocked event of its reason where the label ended:
# "blocked continuations are always reported"
_BLOCKED_ENDS = {"edge_reuse", "edge_reuse_rc", "rejoined_seed"}
# the defaults (§5) of the knobs a structural end needs
_STRUCTURAL_KNOBS = {("branching", "max_label_branches"): 0, ("branching", "max_splits_per_path"): 64,
                     ("branching", "min_successor_labels"): 1, ("branching", "min_successor_fraction"): 0.0,
                     ("bounds", "min_live_labels"): 1}
# the same knobs under the exhaustive preset (§6.9: unlimited branches and splits, no quorum) and in annotate mode
# (§5: no lineage to count, the quorum knobs refused)
_PRESET_KNOBS = {("branching", "max_label_branches"): "unlimited", ("branching", "max_splits_per_path"): "unlimited",
                 ("branching", "min_successor_labels"): 1, ("branching", "min_successor_fraction"): 0.0,
                 ("bounds", "min_live_labels"): 1}


def _requested_strategies(sc: dict, arm: Optional[str]) -> List[Tuple[str, dict]]:
    """The strategies that declare the walk an arm asked for, in order, each (whose, strategy): on the MCP arm the
    overrides of the scenario's MCP arguments (the knobs it sends beside its preset), then the scenario's direct
    request, the walk the preset mirrors; on a direct arm the request."""
    out: List[Tuple[str, dict]] = []   # ("overrides" | "request", strategy)
    if arm == "mcp":
        args = (sc.get("mcp") or {}).get("arguments") if isinstance(sc.get("mcp"), dict) else None
        ov = args.get("overrides") if isinstance(args, dict) else None
        if isinstance(ov, dict):
            out.append(("overrides", ov))
    st = (sc.get("request") or {}).get("strategy") if isinstance(sc.get("request"), dict) else None
    if isinstance(st, dict):
        out.append(("request", st))
    return out


def _declared_knob(sc: dict, group: str, name: str, arm: Optional[str]) -> Tuple[Optional[str], Any]:
    """(whose, value) of strategy.<group>.<name> as the arm's requests declare it (_requested_strategies, whose
    "overrides" or "request"): a knob
    one states, else the default of the first one's exhaustive preset or annotate mode (§6.9); (None, None) where
    none declares it."""
    reqs = _requested_strategies(sc, arm)
    for whose, st in reqs:
        g = st.get(group)
        if isinstance(g, dict) and name in g:
            return whose, g[name]
    for whose, st in reqs:
        if st.get("exhaustive") is True or (st.get("labels") or {}).get("mode") == "annotate":
            return whose, _PRESET_KNOBS.get((group, name))
    return None, None


# The strategy.clamped fields the engine states, each with the knob of the strategy it lowers (src/cli/traverse.cpp,
# the clamps of a /traverse request): bounds.time_budget_ms above the server's max_time_ms (or a budget <= 0 raised
# to it where a seed derives its labels: that budget also bounds the derivation), labels.max_seed_labels above its
# cap, bounds.max_memory_mb and bounds.max_work_units above the server's maxima (or omitted and set to them, stated
# with `requested` "unlimited"). A structural knob is never clamped: out of range it is refused (limit_or_unlimited,
# Strict ranges), so a clamp of one is no engine statement and never stands for the knob.
_CLAMPABLE = {"bounds.time_budget_ms": ("bounds", "time_budget_ms"),
              "labels.max_seed_labels": ("labels", "max_seed_labels"),
              "bounds.max_memory_mb": ("bounds", "max_memory_mb"),
              "bounds.max_work_units": ("bounds", "max_work_units")}
_CLAMPED_FROM_UNLIMITED = {"bounds.max_memory_mb", "bounds.max_work_units"}


def _is_num(v: Any) -> bool:
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def _seeds_derive(sc: dict, arm: Optional[str]) -> bool:
    """Whether a seed of the arm's request derives its labels (none named: SPEC-labeled-traversal-core §5), True
    where the seeds are not known (the MCP arm's arguments, else the scenario's request)."""
    seeds = None
    if arm == "mcp":
        args = (sc.get("mcp") or {}).get("arguments") if isinstance(sc.get("mcp"), dict) else None
        seeds = args.get("seeds") if isinstance(args, dict) else None
    if not isinstance(seeds, list):
        seeds = (sc.get("request") or {}).get("seeds") if isinstance(sc.get("request"), dict) else None
    if not isinstance(seeds, list) or not seeds:
        return True
    return any(not (s.get("labels") if isinstance(s, dict) else None) for s in seeds)


def _scenario_clamp_cap(sc: dict, field: str, where: str = "strategy") -> Optional[float]:
    """The cap of a clampable field the scenario knows: the effective value its own `clamped` expectation of that
    field states (|where|: "strategy" for /traverse's strategy.clamped, "limits" for /resolve's limits.clamped,
    check_clamped's default), the builder's from the host's capabilities snapshot; for /resolve's time budget also
    budget_respected's cap_ms (resolve.time_budget.max_time_ms); else None."""
    for exp in sc.get("expect") or []:
        args = exp.get("args") if isinstance(exp, dict) else None
        if not isinstance(args, dict):
            continue
        if exp.get("check") == "clamped" and args.get("where", "limits") == where and \
                args.get("field") == field and _is_num(args.get("effective")):
            return float(args["effective"])
        if exp.get("check") == "budget_respected" and where == "limits" and field == "bounds.time_budget_ms" and \
                _is_num(args.get("cap_ms")):
            return float(args["cap_ms"])
    return None


def _clamp_problem(sc: dict, ans: dict, c: Any, arm: Optional[str]) -> Optional[str]:
    """Why one strategy.clamped entry cannot be the engine's, or None. The engine states a clamp only of the four
    _CLAMPABLE fields, with `requested` (the request's value, "unlimited" for an omitted budget) and `effective`
    (the server's positive cap), and only as a lowering (the two exceptions: a time budget <= 0 raised to the cap
    where a seed derives its labels, a budget from "unlimited"). Where it is known the entry also agrees with the
    request (a direct arm's request is what the runner sent; the MCP arm's is the service's, so not compared), the
    cap (the scenario's own `clamped` expectation of the field) and the strategy echo, which states the effective
    value the walk ran under."""
    if not isinstance(c, dict):
        return "strategy.clamped entry %r is not an object" % (c,)
    f = c.get("field")
    if f not in _CLAMPABLE:
        return "strategy.clamped names %r, a field the engine never clamps (only %s): no clamp stands for the " \
               "request's value" % (f, ", ".join(sorted(_CLAMPABLE)))
    if "requested" not in c:
        return "strategy.clamped %s states no requested value" % f
    req, eff = c["requested"], c.get("effective")
    if not _is_num(eff) or eff <= 0:
        return "strategy.clamped %s: effective %r is no server cap (a positive number)" % (f, eff)
    if req == "unlimited":
        if f not in _CLAMPED_FROM_UNLIMITED:
            return "strategy.clamped %s: requested \"unlimited\", which only an omitted budget states" % f
    elif not _is_num(req):
        return "strategy.clamped %s: requested %r is not a number" % (f, req)
    elif not (eff < req or (f == "bounds.time_budget_ms" and req <= 0 and _seeds_derive(sc, arm))):
        return "strategy.clamped %s: %r -> %r is no lowering to the server's cap%s" % (
            f, req, eff, " (a budget <= 0 is raised only where a seed derives its labels)"
            if f == "bounds.time_budget_ms" and req <= 0 else "")
    group, name = _CLAMPABLE[f]
    if arm != "mcp":
        st = (sc.get("request") or {}).get("strategy") if isinstance(sc.get("request"), dict) else None
        g = st.get(group) if isinstance(st, dict) else None
        if isinstance(g, dict) and name in g:
            if g[name] != req:
                return "strategy.clamped %s: requested %r, the request sent %r" % (f, req, g[name])
        elif f in _CLAMPED_FROM_UNLIMITED and req != "unlimited":
            return "strategy.clamped %s: requested %r, the request omitted the budget (\"unlimited\")" % (f, req)
    cap = _scenario_clamp_cap(sc, f)
    if cap is not None and eff != cap:
        return "strategy.clamped %s: effective %r, the server's cap is %r" % (f, eff, cap)
    body = ans.get("body") if isinstance(ans.get("body"), dict) else {}
    echo = body.get("strategy") if isinstance(body.get("strategy"), dict) else {}
    eg = echo.get(group)
    if isinstance(eg, dict) and _is_num(eg.get(name)) and eg[name] != eff:
        return "strategy.clamped %s: effective %r, the strategy echo ran under %r" % (f, eff, eg[name])
    return None


def _clamp_problems(sc: dict, ans: dict, arm: Optional[str]) -> List[str]:
    """Why the answer's strategy.clamped cannot be the engine's (_clamp_problem per entry; a list, one entry per
    field), and, on a direct arm, a clampable knob the request sent that the strategy echo states otherwise without a
    clamp naming it: the engine echoes the request's value unless it clamped it (traverse.cpp, the clamps of a
    request), and the walk's checks read the echo (a time budget lowered unstated would excuse a premature time
    stop, one raised unstated an overrun). The MCP arm's request is the service's: not compared."""
    body = ans.get("body") if isinstance(ans.get("body"), dict) else {}
    echo = body.get("strategy") if isinstance(body.get("strategy"), dict) else {}
    lst = echo.get("clamped")
    if lst is not None and not isinstance(lst, list):
        return ["strategy.clamped %r is not a list" % (lst,)]
    out, seen = [], set()
    named = {c.get("field") for c in lst or [] if isinstance(c, dict)}   # stated, a problem reported if invalid
    for c in lst or []:
        bad = _clamp_problem(sc, ans, c, arm)
        if bad:
            out.append(bad)
            continue
        if c["field"] in seen:
            out.append("strategy.clamped states %s twice" % c["field"])
        seen.add(c["field"])
    st = (sc.get("request") or {}).get("strategy") if isinstance(sc.get("request"), dict) else None
    if arm != "mcp" and isinstance(st, dict):
        for f, (group, name) in sorted(_CLAMPABLE.items()):
            g, eg = st.get(group), echo.get(group)
            if f in named or not isinstance(g, dict) or not isinstance(eg, dict):
                continue
            sent, ran = g.get(name), eg.get(name)
            if _is_num(sent) and _is_num(ran) and ran != sent:
                out.append("strategy.%s %r, the request sent %r, without a clamp of it" % (f, ran, sent))
    return out


def _resolve_clamp_problems(sc: dict, body: dict, arm: Optional[str]) -> List[str]:
    """Why a /resolve answer's limits.clamped cannot be the engine's. The route clamps one field, a time budget above
    the server's resolve.time_budget.max_time_ms, lowered to it (traverse.cpp process_resolve_request: `requested` the
    value as given, `effective` the cap; a budget not above the finalisation reserve is refused, never raised), and
    states limits.time_budget_ms, the budget the work ran under. So an entry names bounds.time_budget_ms with a
    numeric `requested` above its positive `effective`, once; where it is known it agrees with the request (a direct
    arm's; the MCP arm's request is the service's, not compared), the cap (the scenario's `clamped` expectation or
    budget_respected's cap_ms) and limits.time_budget_ms. A direct arm's limits.time_budget_ms below the request's
    budget without that clamp is a lowering the engine never makes unstated."""
    lim = body.get("limits") if isinstance(body, dict) else None
    if not isinstance(lim, dict):
        return []
    field = "bounds.time_budget_ms"
    req = sc.get("request") if isinstance(sc.get("request"), dict) else {}
    sent = (req.get("bounds") or {}).get("time_budget_ms") if isinstance(req.get("bounds"), dict) else None
    ran = lim.get("time_budget_ms")
    lst = lim.get("clamped")
    if lst is not None and not isinstance(lst, list):
        return ["limits.clamped %r is not a list" % (lst,)]
    out, seen = [], set()
    for c in lst or []:
        if not isinstance(c, dict):
            out.append("limits.clamped entry %r is not an object" % (c,))
            continue
        f = c.get("field")
        if f != field:
            out.append("limits.clamped names %r, a field /resolve never clamps (only %s)" % (f, field))
            continue
        if f in seen:
            out.append("limits.clamped states %s twice" % f)
            continue
        seen.add(f)
        if "requested" not in c:
            out.append("limits.clamped %s states no requested value" % f)
            continue
        rq, eff = c["requested"], c.get("effective")
        if not _is_num(eff) or eff <= 0:
            out.append("limits.clamped %s: effective %r is no server cap (a positive number)" % (f, eff))
        elif not _is_num(rq):
            out.append("limits.clamped %s: requested %r is not a number" % (f, rq))
        elif not eff < rq:
            out.append("limits.clamped %s: %r -> %r is no lowering to the server's cap (/resolve refuses a budget "
                       "it cannot keep, it never raises one)" % (f, rq, eff))
        elif arm != "mcp" and rq != sent:
            out.append("limits.clamped %s: requested %r, the request sent %r" % (f, rq, sent))
        elif _scenario_clamp_cap(sc, f, "limits") not in (None, eff):
            out.append("limits.clamped %s: effective %r, the server's cap is %r" % (
                f, eff, _scenario_clamp_cap(sc, f, "limits")))
        elif _is_num(ran) and ran != eff:
            out.append("limits.clamped %s: effective %r, the work ran under limits.time_budget_ms %r" % (f, eff, ran))
    if arm != "mcp" and field not in seen and _is_num(sent) and _is_num(ran) and ran != sent:
        out.append("limits.time_budget_ms %r, the request sent %r, without a clamp of it" % (ran, sent))
    return out


def _clamped_knob(sc: dict, ans: dict, group: str, name: str, arm: Optional[str] = None) -> Tuple[bool, Any]:
    """(True, the effective value) where the answer's single strategy.clamped entry of strategy.<group>.<name> is
    one the engine can state (_clamp_problem: a _CLAMPABLE field, lowered to its cap); else (False, None): a clamp
    of another knob, or one that cannot have happened, never stands for the request's value."""
    field = "%s.%s" % (group, name)
    if field not in _CLAMPABLE:
        return False, None
    body = ans.get("body") if isinstance(ans.get("body"), dict) else {}
    echo = body.get("strategy") if isinstance(body.get("strategy"), dict) else {}
    lst = echo.get("clamped")
    mine = [c for c in lst if isinstance(c, dict) and c.get("field") == field] if isinstance(lst, list) else []
    if len(mine) != 1 or _clamp_problem(sc, ans, mine[0], arm):
        return False, None
    return True, mine[0]["effective"]


def _strategy_knob(sc: dict, ans: dict, group: str, name: str, arm: Optional[str] = None) -> Any:
    """The value of strategy.<group>.<name> the walk is held to (SPEC-labeled-traversal-core §5): a clamp's
    effective value (strategy.clamped, only a clamp the engine can state: _clamped_knob, so never one of a
    structural knob), else what the arm's requests declare (_declared_knob: on the MCP arm its overrides, then the
    scenario's direct request, whose walk the preset mirrors), else the answer's echo, else the default (the echo's
    preset or mode, else §5's). The declaration comes before the echo: a knob a service dropped or rewrote, a
    preset's included, is held to what the scenario asked, never excused by its own echo or a clamp it states
    (_preset_mismatch names such a difference, _clamp_problems such a clamp)."""
    clamped, value = _clamped_knob(sc, ans, group, name, arm)
    if clamped:
        return value
    whose, value = _declared_knob(sc, group, name, arm)
    if whose is not None:
        return value
    body = ans.get("body") if isinstance(ans.get("body"), dict) else {}
    echo = body.get("strategy") if isinstance(body.get("strategy"), dict) else {}
    g = echo.get(group)
    if isinstance(g, dict) and name in g:
        return g[name]
    if echo.get("exhaustive") is True or (echo.get("labels") or {}).get("mode") == "annotate":
        return _PRESET_KNOBS.get((group, name))
    return _STRUCTURAL_KNOBS.get((group, name))


def _preset_mismatch(sc: dict, ans: dict, arm: Optional[str]) -> List[str]:
    """On the MCP arm, every structural walk knob whose echo differs from what the scenario declares for it (its MCP
    overrides, else its direct request, _declared_knob) without a clamp: the service walked another strategy than
    the scenario's, a preset or override it changed or dropped (a deliberate preset change needs the manifest to
    follow it)."""
    if arm != "mcp":
        return []
    body = ans.get("body") if isinstance(ans.get("body"), dict) else {}
    echo = body.get("strategy") if isinstance(body.get("strategy"), dict) else {}
    out = []
    for group, name in _STRUCTURAL_KNOBS:
        g = echo.get(group)
        if not isinstance(g, dict) or name not in g or _clamped_knob(sc, ans, group, name, arm)[0]:
            continue
        whose, declared = _declared_knob(sc, group, name, arm)
        if whose == "overrides" and g[name] != declared:
            out.append("the service walked with %s.%s %r, the scenario's MCP overrides state %r: an override not "
                       "applied" % (group, name, g[name], declared))
        elif whose == "request" and g[name] != declared:
            out.append("the service walked with %s.%s %r, the scenario declares %r: a preset mismatch" % (
                group, name, g[name], declared))
    return out


def _route_splits(armo: dict) -> Callable[[Any], Tuple[int, int]]:
    """segment id -> (the most ambiguous splits, the most splits) any root-to-segment route passes before reaching
    the segment, through any parent at a merge. A split sits at the end of its segment (SPEC §7.1 splits[].segment;
    a graphlet's G split flag); a split whose kind is not stated counts as ambiguous."""
    segs = {s.get("id"): s for s in armo.get("segments") or [] if isinstance(s, dict)}
    amb: Dict[Any, bool] = {}
    for sp in armo.get("splits") or []:
        if not isinstance(sp, dict):
            continue
        kind = sp.get("kind")
        amb[sp.get("segment")] = kind != "divergence" if kind is not None else sp.get("ambiguous") is not False
    memo: Dict[Any, Tuple[int, int]] = {}

    def most(sid: Any) -> Tuple[int, int]:
        if sid in memo:
            return memo[sid]
        memo[sid] = (0, 0)   # a cycle (never in a valid answer) counts nothing
        a = n = 0
        for par in (segs.get(sid) or {}).get("parents") or []:
            pa, pn = most(par)
            if par in amb:
                pa, pn = pa + (1 if amb[par] else 0), pn + 1
            a, n = max(a, pa), max(n, pn)
        memo[sid] = (a, n)
        return memo[sid]
    return most


def _finite(v: Any) -> bool:
    return isinstance(v, int) and not isinstance(v, bool)


def _arm_index(armo: dict) -> dict:
    """What _structural_end_problems reads of an arm, built once per arm (every path of a wide trie would otherwise
    rebuild it: quadratic in the paths): the segments by id, the route split counts (_route_splits) and the runs'
    branches by (segment, to_bp, label)."""
    runs: Dict[Tuple[Any, Any], Dict[Any, Any]] = {}
    for x in armo.get("runs") or []:
        if isinstance(x, dict):
            runs.setdefault((x.get("segment"), x.get("to_bp")), {})[x.get("label")] = x.get("branches")
    return {"segments": {x.get("id"): x for x in armo.get("segments") or [] if isinstance(x, dict)},
            "splits": _route_splits(armo), "runs": runs}


def _structural_end_problems(sc: dict, ans: dict, armo: dict, p: dict, index: Optional[dict] = None,
                             arm: Optional[str] = None) -> List[str]:
    """Why the structural ends stated at path |p|'s leaf cannot have happened, judged from the request, the strategy
    (_strategy_knob) and the arm (SPEC-labeled-traversal-core §6.4, §6.6, §6.7):
    - a structural reason the path counts (end_reasons) has a label_end event of it on its leaf at the path's end
      (judged where the leaf states its events);
    - `branch`: the lineage was ambiguous with its allowance spent (branches + 1 > max_label_branches, and a lineage
      never has more): max_label_branches finite, the label's branches (end_labels, or runs[] at the leaf) equal to
      it, and two or more structural successors where the event states them; without a quorum knob also a route to
      the leaf passing at least that many ambiguous splits (each branch taken is a split there). Under a quorum
      knob a branch can be taken without a split: the walker raises the lineage's branches on every admissible
      successor before the quorum drops the minority ones, and one followed successor records no split
      (walker.cpp, the ambiguity before the quorum);
    - `minority`: a quorum knob set (min_successor_labels >= 2 or min_successor_fraction > 0); `below_min_labels`:
      min_live_labels > 1; `split_limit`: max_splits_per_path finite and a route to the leaf through that many
      splits;
    - `edge_reuse`, `edge_reuse_rc`, `rejoined_seed`: a blocked event of that reason on the leaf at the path's end,
      also where the path states the re-entry only as its path_reason (annotate mode: no lineage, end_reasons
      empty); any structural end, counted or a path_reason, needs the leaf's events list (Missing without it).
    |index| is the arm's _arm_index, built here when not given; |arm| names the arm whose request holds the knobs."""
    out: List[str] = []
    chain = p.get("segments") or []
    if not chain:
        return out
    index = index if index is not None else _arm_index(armo)
    leaf = index["segments"].get(chain[-1])
    if leaf is None:
        return out
    L = p.get("length_bp")
    where = "path %s at %s bp" % (p.get("id"), L)
    claimed = set(p.get("end_reasons") or {}) & _STATED_STRUCTURAL_ENDS
    # annotate mode ends a path at a re-entry by its path_reason alone: no lineage, end_reasons empty (walker.cpp,
    # the structural end of an annotate walk: finish_path with the block's reason, its blocked event pushed there)
    by_reason = {p.get("path_reason")} & _STATED_STRUCTURAL_ENDS
    if (claimed or by_reason) and not isinstance(leaf.get("events"), list):
        # the evidence of a structural end is its leaf's events (§7.0: below complete_to_bp every decision is
        # reported; traverse.cpp writes `events` on every segment of a tree or full detail, a graphlet states them in
        # its E and R records, graphlet_arms lists them on every segment): without the list the end is not shown,
        # never excused (fail on a direct arm, not exposed on the MCP's)
        raise Missing("segments[%s].events of %s, a %s end (%s)" % (
            chain[-1], where, "/".join(sorted(claimed)) if claimed else "path_reason %s" % p.get("path_reason"),
            "absent" if "events" not in leaf else "%s, not a list" % type(leaf["events"]).__name__))
    evs = [x for x in _seg_events(leaf) if x.get("at_bp") == L]
    ends = [x for x in evs if x.get("type") == "label_end"]

    def named(ev: dict) -> Any:
        text = ev.get("text")
        return text if ev.get("reason") == "branch" and text in _QUORUM_ENDS else ev.get("reason")
    if isinstance(leaf.get("events"), list):
        stated = {named(x) for x in ends}
        stated |= {"branch" for x in stated if x in _QUORUM_ENDS}
        for reason in sorted(claimed - stated):
            out.append("%s ends by %s without a label_end event of it on its leaf" % (where, reason))
        # a re-entry path_reason states its blocked event where the path ends (§6.6: "blocked continuations are
        # always reported"), in any mode
        for reason in sorted(by_reason & _BLOCKED_ENDS):
            if not any(x.get("type") == "blocked" and x.get("reason") == reason for x in evs):
                out.append("%s ends by path_reason %s without a blocked event of it on its leaf" % (where, reason))
    mlb = _strategy_knob(sc, ans, "branching", "max_label_branches", arm)
    msp = _strategy_knob(sc, ans, "branching", "max_splits_per_path", arm)
    msl = _strategy_knob(sc, ans, "branching", "min_successor_labels", arm)
    msf = _strategy_knob(sc, ans, "branching", "min_successor_fraction", arm)
    mll = _strategy_knob(sc, ans, "bounds", "min_live_labels", arm)
    quorum = (_finite(msl) and msl >= 2) or (isinstance(msf, (int, float)) and msf > 0) or (_finite(mll) and mll > 1)
    amb_n, split_n = index["splits"](chain[-1])
    lineage = dict(index["runs"].get((chain[-1], L)) or {})
    lineage.update({x.get("label"): x.get("branches") for x in p.get("end_labels") or [] if isinstance(x, dict)})
    for ev in ends:
        what = named(ev)
        lab = "%s label %s" % (where, ev.get("label"))
        if what == "branch":
            br = lineage.get(ev.get("label"))
            if not _finite(mlb):
                out.append("%s: a branch end with max_label_branches %r" % (lab, mlb))
            elif _finite(br) and br != mlb:
                out.append("%s: a branch end after %d branch(es) of max_label_branches %d" % (lab, br, mlb))
            elif amb_n < mlb and not quorum:
                out.append("%s: a branch end after %d ambiguous split(s) on any route to it, max_label_branches %d"
                           % (lab, amb_n, mlb))
            elif _finite(ev.get("structural_successors")) and ev["structural_successors"] < 2:
                out.append("%s: a branch end at %d structural successor(s)" % (lab, ev["structural_successors"]))
        elif what == "minority":
            if not ((_finite(msl) and msl >= 2) or (isinstance(msf, (int, float)) and msf > 0)):
                out.append("%s: a minority end with min_successor_labels %r and min_successor_fraction %r" % (
                    lab, msl, msf))
        elif what == "below_min_labels":
            if not (_finite(mll) and mll > 1):
                out.append("%s: a below_min_labels end with min_live_labels %r" % (lab, mll))
        elif what == "split_limit":
            if not _finite(msp) or split_n < msp:
                out.append("%s: a split_limit end after %d split(s) on any route to it, max_splits_per_path %r" % (
                    lab, split_n, msp))
        elif what in _BLOCKED_ENDS:
            if not any(x.get("type") == "blocked" and x.get("reason") == what for x in evs):
                out.append("%s: ended by %s without a blocked event of it there" % (lab, what))
    return out


def _seed_frame(t: dict, c: dict, side: str, L: int) -> Tuple[Callable[[int], str], int, int]:
    """For one copy |c| of the truth window: (continuation(n), room in the flank, distance to the record's end) of
    the arm |side| of a seed that is the window ('+') or its reverse complement ('-')."""
    s0 = c["start"] - c["flank_start"]
    fl = c["flank"]
    rec_len, start = t["record_length"], c["start"]
    minus = t.get("strand", "+") == "-"
    walks_right_on_record = (side == "right") != minus
    # a '-' seed's right arm reads the record leftwards on the reverse strand: rc(w + Y) = rc(Y) + rc(w)
    if walks_right_on_record:
        def cont(n: int) -> str:
            x = fl[s0 + L:s0 + L + n]
            return revcomp(x) if minus else x
        room, end = len(fl) - (s0 + L), rec_len - (start + L)
    else:
        def cont(n: int) -> str:
            x = fl[max(0, s0 - n):s0]
            return revcomp(x) if minus else x
        room, end = s0, start
    return cont, room, end


def _regime(ans: dict, r: dict) -> str:
    """The graph regime of the answer's index (basic, canonical, primary): the body's capabilities, else the graphlet's
    H record; basic when neither states it."""
    body = ans.get("body") if isinstance(ans.get("body"), dict) else {}
    reg = (body.get("capabilities") or {}).get("regime") if isinstance(body.get("capabilities"), dict) else None
    if not reg and isinstance(r.get("graphlet"), str):
        head = r["graphlet"].split("\n", 1)[0].split(" ")
        reg = head[4] if len(head) > 4 and head[0] == "H" else None
    return str(reg or "basic")


def _first_block(seed: str, out: str, side: str, canonical: bool, k: int = K,
                 hairpins: str = "skip") -> Optional[Tuple[int, str]]:
    """(d, reason) of the first outward depth d at which the step from the seed's arm along |out| (the record's bases
    in walking order) to out[d] is inadmissible by SPEC-labeled-traversal-core §6.5 / §6.6, or None: in
    canonical/primary a hairpin step under |hairpins| skip (its (k+1)-mer its own reverse complement, at even k also
    a step into or out of a palindromic k-mer), first; else rejoined_seed (the step's target k-mer is a seed k-mer,
    either orientation in canonical/primary), else edge_reuse_rc / edge_reuse (the step's (k+1)-mer, or in
    canonical/primary its reverse complement when that is another (k+1)-mer, is an earlier step of the walk: a
    palindromic one reused is edge_reuse, as the engine keys it); precedence rejoined_seed > edge_reuse_rc >
    edge_reuse. The oracle's Walker.block decides each step alike (tests/test_oracle_traverse.py FirstBlockParity)."""
    if len(seed) < k:
        return None
    seedk = set(seed[j:j + k] for j in range(len(seed) - k + 1))
    if canonical:
        seedk |= set(revcomp(x) for x in seedk)
    used: set = set()
    win = seed[-k:] if side == "right" else seed[:k]   # the walk's node, natural orientation
    for d, b in enumerate(out):
        kmer, edge = (win[1:] + b, win + b) if side == "right" else (b + win[:-1], b + win)
        rc = revcomp(edge)
        if canonical and hairpins == "skip" and (
                edge == rc or (k % 2 == 0 and (win == revcomp(win) or kmer == revcomp(kmer)))):
            return d, "hairpin"
        if kmer in seedk:
            return d, "rejoined_seed"
        if canonical and rc != edge and rc in used:
            return d, "edge_reuse_rc"
        if edge in used:
            return d, "edge_reuse"
        used.add(edge)
        win = kmer
    return None


def _record_view(t: dict, c: dict, L: int) -> Tuple[str, int, int, Callable[[int], int]]:
    """For one copy |c| of the truth window: (the truth flank in the seed's orientation O, the seed's span [p0, p1) in
    O, the record coordinate of O[x])."""
    fl = c["flank"]
    s0 = c["start"] - c["flank_start"]
    if t.get("strand", "+") == "-":
        return revcomp(fl), len(fl) - (s0 + L), len(fl) - s0, lambda x: c["flank_start"] + len(fl) - 1 - x
    return fl, s0, s0 + L, lambda x: c["flank_start"] + x


def _occurrences(text: str, word: str) -> List[int]:
    out, j = [], text.find(word)
    while j >= 0:
        out.append(j)
        j = text.find(word, j + 1)
    return out


_MAX_REPEAT_EXITS = 64


# the path end of a computed block (_computed_block): a skipped hairpin ends the path as dead_end (SPEC §13), every
# other block under its own name
_GAP_ENDS = {"hairpin": "dead_end"}


def _gap_stated(gap: dict, ends: set, armo: dict, routes: List[Tuple[dict, str]], best: int) -> bool:
    """Whether the answer states the computed block |gap| where the longest routes following the record end (|ends|
    their leaves' end reasons, _leaf_ends): its end reason (_GAP_ENDS), and for a hairpin the hairpin itself, the
    label's `hairpin` end (SPEC-labeled-traversal-core §13 note 2: dead_end with text hairpin; a graphlet's Dh) or
    the skipped `hairpin` event at the gap's depth with the record's base on such a route's leaf (§6.5: "a hairpin
    event is recorded"; annotate mode states no label end). A bare dead_end there is a walk that stopped short (L3 of
    the verification of SD-025)."""
    if _GAP_ENDS.get(gap["reason"], gap["reason"]) not in ends:
        return False
    if gap["reason"] != "hairpin" or "hairpin" in ends:
        return True
    segs = dict((x.get("id"), x) for x in armo.get("segments") or [])
    for p, f in routes:
        if len(f) != best or not p.get("segments"):
            continue
        for ev in _seg_events(segs.get(p["segments"][-1])):
            if ev.get("type") == "hairpin" and ev.get("followed") is not True and ev.get("at_bp") == gap["d"] \
                    and ev.get("char") == gap["char"]:
                return True
    return False


def _computed_block(sc: dict, ans: dict, r: dict, i: int, t: dict, c: dict, side: str, L: int,
                    arm: Optional[str] = None) -> Optional[dict]:
    """The first depth at which the record's own continuation from the seed is inadmissible (_first_block), computed
    from the truth's flank and the request's seed, never read from the answer, under the hairpin policy the arm's
    requests declare (_strategy_knob with |arm|: on the MCP arm its overrides): {d, reason, char (the record's base
    there), exits: [(the record coordinate of the exit base, the outward bases from it)] (the junction's other
    successors in the flank: where the walk leaves a collapsed repeat), reaches_end (the flank holds the record's end
    in the walking direction)}, or None."""
    seeds = sc["request"].get("seeds") or []
    seed = str((seeds[i] if i < len(seeds) and isinstance(seeds[i], dict) else {}).get("sequence") or "").upper()
    O, p0, p1, coord = _record_view(t, c, L)
    out = O[p1:] if side == "right" else O[:p0][::-1]
    blk = _first_block(seed, out, side, _regime(ans, r) in ("canonical", "primary"),
                       hairpins=_strategy_knob(sc, ans, "branching", "hairpins", arm) or "skip")
    if blk is None or blk[0] >= len(out):
        return None
    d, why = blk
    own = out[d]
    exits = []
    if side == "right":
        junction = (seed + out[:d])[-(K - 1):]
        for x in _occurrences(O, junction):
            if x + K - 1 < len(O) and O[x + K - 1] != own:
                exits.append((coord(x + K - 1), O[x + K - 1:]))
        x_end = len(O) - 1
    else:
        junction = (out[:d][::-1] + seed)[:K - 1]
        for x in _occurrences(O, junction):
            if x >= 1 and O[x - 1] != own:
                exits.append((coord(x - 1), O[:x][::-1]))
        x_end = 0
    return {"d": d, "reason": why, "char": own, "out": out, "exits": exits[:_MAX_REPEAT_EXITS],
            "reaches_end": coord(x_end) in (0, t["record_length"] - 1)}


def _blocked_walk(sc: dict, ans: dict, r: dict, i: int, t: dict, armo: dict, side: str, copies: List[dict], L: int,
                  labels: Optional[set], exact: Optional[dict], record_scope: bool,
                  arm: Optional[str] = None) -> Optional[Tuple[str, str]]:
    """G28 (VERIFY-2 item 2, the strict rule): an arm with no route that follows the record for its whole length
    still follows it when the record's own continuation is blocked (SPEC-labeled-traversal-core §6.6: a step into a
    seed node or along a used edge is not taken but reported as a `blocked` event): some route (any parent at a
    merge; the label rule as _record_routes) spells the record for exactly its first d bases, where d is the first
    inadmissible depth computed from the truth's flank and the request's seed (_computed_block); one of its segments
    states the blocked event at at_bp d with the record's base and the computed reason; and past d it spells the
    record from a collapsed repeat's exit (an occurrence of the junction in the flank followed by another base) for
    every base the flank holds. (outcome, detail): pass; inconclusive when d lies before the truth's min_bp (a truth
    gap), when such a route runs past the flank (not verifiable), or when the search stopped at its bound; fail when
    it runs past the record's end with the record as its only label; None when no copy has such a route. A block
    whose repeat has no exit inside the truth's flank (A2-04 of the review of lane S5b): a route that spells the
    record's first d bases, states the block and continues is INCONCLUSIVE (the continuation cannot be verified),
    one that ends at the block passes, none: None."""
    unverifiable = None   # a copy whose repeat exit lies outside the flank: the answer when no other copy decides
    for c in copies:
        blk = _computed_block(sc, ans, r, i, t, c, side, L, arm)
        if blk is None:
            continue
        if not blk["exits"]:
            found = _blocked_without_exit(sc, armo, side, blk, labels, exact)
            if found is not None and found[0] != INCONCLUSIVE:
                return found
            unverifiable = unverifiable or found
            continue
        d, out = blk["d"], blk["out"]
        frames = []
        for pos, tail in blk["exits"]:
            walk = out[:d] + tail
            frames.append(((lambda n, w=walk: w[:n]) if side == "right" else (lambda n, w=walk: w[:n][::-1]),
                           len(walk), len(walk)))
        routes, beyond, capped = _record_routes(armo, side, frames, labels=labels, event=(d, blk["char"], blk["reason"]))
        what = "the record's own continuation blocked at %d bp by %s (computed from the record and the seed)" % (
            d, blk["reason"])
        if capped:
            return INCONCLUSIVE, "%s; the route search stopped at %d states" % (what, _MAX_ROUTE_STATES)
        if routes:
            longest = max(len(f) for _, f in routes)
            exit_at = sorted(set(pos for pos, tail in blk["exits"] if any(
                f == (out[:d] + tail)[:len(f)] if side == "right" else f[::-1] == (out[:d] + tail)[:len(f)]
                for _, f in routes)))
            if exact is not None and d < exact.get("min_bp", 0):
                return INCONCLUSIVE, "%s, before the truth's minimum of %d bp: a truth gap" % (what, exact["min_bp"])
            return PASS, "followed for %d bp, then blocked by %s (computed from the record and the seed, stated at " \
                         "%d bp); %d route(s) continue along the record from its exit at %s (longest %d bp)" % (
                             d, blk["reason"], d, len(routes), exit_at[:3], longest)
        if beyond:
            if record_scope and blk["reaches_end"]:
                return FAIL, "%s; %d route(s) continue from a repeat exit past the record's end" % (what, beyond)
            return INCONCLUSIVE, "%s; %d route(s) continue from a repeat exit past the truth's flank (not " \
                                 "verifiable)" % (what, beyond)
    return unverifiable


def _blocked_without_exit(sc: dict, armo: dict, side: str, blk: dict, labels: Optional[set],
                          exact: Optional[dict]) -> Optional[Tuple[str, str]]:
    """_blocked_walk for a computed block at d whose collapsed repeat has no exit inside the truth's flank (A2-04):
    the routes that spell the record's first d bases and state the blocked event at d (any parent at a merge, the
    label rule); one that continues past d leaves from an exit the flank does not hold: INCONCLUSIVE (not
    verifiable); one that ends at d: PASS; a block before the truth's min_bp: INCONCLUSIVE (a truth gap); none:
    None (the arm does not follow the record)."""
    d, out = blk["d"], blk["out"]
    walk = out[:d]
    frames = [((lambda n, w=walk: w[:n]) if side == "right" else (lambda n, w=walk: w[:n][::-1]), d, d)]
    routes, beyond, capped = _record_routes(armo, side, frames, labels=labels, event=(d, blk["char"], blk["reason"]))
    what = "the record's own continuation blocked at %d bp by %s (computed from the record and the seed)" % (
        d, blk["reason"])
    if capped:
        return INCONCLUSIVE, "%s; the route search stopped at %d states" % (what, _MAX_ROUTE_STATES)
    if not routes and not beyond:
        return None
    if exact is not None and d < exact.get("min_bp", 0):
        return INCONCLUSIVE, "%s, before the truth's minimum of %d bp: a truth gap" % (what, exact["min_bp"])
    if beyond:
        return INCONCLUSIVE, "%s, stated at %d bp; %d route(s) continue from the repeat's exit, which is not in the " \
                             "truth's flank: not verifiable" % (what, d, beyond)
    return PASS, "followed for %d bp, then blocked by %s (computed from the record and the seed, stated at %d bp); " \
                 "the route ends at the block" % (d, blk["reason"], d)


def check_walk_in_record(sc, ans, a, arm):
    """The walk follows the record: per requested arm, some root-to-leaf route through the segments, taking any
    parent at a merge (_record_routes: a merged walk displays first parents only), spells the record's own
    continuation from the seed (for its whole length; with any_copy, the continuation of any copy of a repeated seed;
    a '-' truth: the seed is the reverse complement of the record's window), the root segment agrees with it, and the
    longest such route reaches its minimum; a route's end reasons are its leaf's (_leaf_ends). A route search that
    stops at its bound before it shows the minimum is inconclusive. The minimum is truth.walk[side].min_bp where the
    builder computed it exactly from the record (header-label walks on whole records: the first base where the
    record offers a second successor on either strand, the first reused edge, the record's end, max_extension_bp),
    and a shorter route fails unless a stated resource limit stopped the arm (inconclusive). Otherwise it is
    min(max_extension_bp, the record's end, the truth's flank), and a shorter route may also end at a stated branch
    or re-entry (_STATED_STRUCTURAL_ENDS: other records of a column or taxa can branch where this record does not).
    With the record's header as the only label (label_scope "record") no flank may run past the record's end. expect_branching: some arm states a
    split or a branch event, at or before truth.walk[side].branch_by_bp where the builder stated it (a copy's
    divergence or a second successor on the record's forward strand). prefix_only: no minimum (a walk a budget
    cut); if_walked: a seed that was not walked (an allowed failed outcome) is inconclusive here."""
    i = a.get("seed", 0)
    r = _result(ans, i)
    if "arms" not in r:
        if a.get("if_walked"):
            return INCONCLUSIVE, "seed %d not walked (outcome %r): nothing to follow" % (i, (r.get("outcome") or {}).get("walks"))
        return FAIL, "seed %d not walked (outcome %r, error %r)" % (i, r.get("outcome"), r.get("error"))
    rejected = _mgt_rejected(r)
    if rejected:
        return FAIL, rejected
    t = sc["truth"] if "truths" not in sc else sc["truths"][i]
    L = t["length"]
    copies = [{"start": t["start"], "flank": t["flank"], "flank_start": t["flank_start"]}]
    if a.get("any_copy"):
        copies += t.get("copies", [])
    maxbp = ((sc["request"].get("strategy") or {}).get("bounds") or {}).get("max_extension_bp")
    probs, notes, stopped, searched = [], [], [], []
    unshown = []   # what the answer does not hold to judge an arm (Missing, once nothing failed)
    mismatch = _preset_mismatch(sc, ans, arm)
    branching = treeless = False
    for side in a.get("arms") or _arms_of(sc, i):
        armo = r["arms"].get(side)
        if armo is None:
            probs.append("arm %s absent" % side)
            continue
        if "segments" not in armo and "paths" not in armo:
            # an arm without its tree (an MCP answer without its graphlet, a summary detail): the walk is not shown,
            # which says nothing of its paths (not "no path")
            unshown.append("results[%d].arms.%s.segments and paths (no graphlet or tree)" % (i, side))
            treeless = True
            continue
        first_branch = _first_branch_bp(armo)
        if first_branch is not None:
            branching = True
            by = ((t.get("walk") or {}).get(side) or {}).get("branch_by_bp")
            if a.get("expect_branching") and by is not None and first_branch > by:
                probs.append("%s: the first split or branch event at %d bp, the record branches by %d bp" % (
                    side, first_branch, by))
        frames = [_seed_frame(t, c, side, L) for c in copies]
        own_end = frames[0][2]
        exact = (t.get("walk") or {}).get(side)
        if exact is not None and "min_bp" in exact:
            expected = exact["min_bp"]
        else:
            exact = None
            expected = min(x for x in (maxbp, own_end, frames[0][1]) if x is not None)
        if a.get("prefix_only"):
            expected = 0
        trig = (armo.get("cap_trigger") or {}).get("reason")
        fl = _arm_flanks(armo, side)
        if not fl:
            if armo.get("status") == "truncated" and trig in _RESOURCE_ENDS:
                stopped.append("%s: no path, stopped by %s" % (side, trig))
            else:
                probs.append("%s: no path (status %s)" % (side, armo.get("status")))
            continue
        if a.get("label_scope") == "record":
            for p, f in fl:
                if len(f) > own_end:
                    probs.append("%s path %s runs %d bp, past the record's end (%d bp)" % (
                        side, p.get("id"), len(f), own_end))
        # the routes that follow the record, through any parent at a merge, on segments carrying its label where they
        # state labels (_record_routes, _route_labels)
        labels, seq_only = _route_labels(r, t, armo)
        how = " [%s]" % seq_only if seq_only else ""
        following, beyond, capped = _record_routes(armo, side, frames, labels=labels)
        if a.get("label_scope") == "record":
            following = [(p, f) for p, f in following if len(f) <= own_end]
        if not following and not capped:
            blocked = _blocked_walk(sc, ans, r, i, t, armo, side, copies, L, labels, exact,
                                    a.get("label_scope") == "record", arm)
            if blocked is not None:
                kind, text = blocked
                (notes if kind == PASS else searched if kind == INCONCLUSIVE else probs).append(
                    "%s: %s%s" % (side, text, how))
                continue
        if not following:
            lens = [len(f) for _, f in fl]
            if capped:
                searched.append("%s: the route search stopped at %d visits of the segment graph before a route "
                                "following the record was found (paths of %s bp)" % (side, _MAX_ROUTE_STATES, lens[:6]))
                continue
            probs.append("%s: no route of the arm's %d segment(s) follows the record, any parent taken at a merge "
                         "(%d path(s) of %s bp%s)%s" % (side, len(armo.get("segments") or []), len(fl), lens[:6],
                                                        "; %d route(s) beyond the truth flank" % beyond if beyond else "",
                                                        how))
            continue
        for sg in [x for x in armo.get("segments", []) if not x.get("parents")]:
            rs = sg.get("sequence", "")
            if not any(room >= len(rs) and rs == cont(len(rs)) for cont, room, _ in frames):
                probs.append("%s: the root segment disagrees with the record" % side)
        best = max(len(f) for _, f in following)
        shown = set(f for _, f in fl)
        merged = "" if any(f in shown for _, f in following if len(f) == best) else \
            ", through a non-first parent of a merge"
        if best >= expected:
            notes.append("%s: %d route(s) follow the record (longest %d bp of %d expected%s)%s" % (
                side, len(following), best, expected, merged, how))
            continue
        if capped:
            searched.append("%s: the record followed for %d of %d bp when the route search stopped at %d visits" % (
                side, best, expected, _MAX_ROUTE_STATES))
            continue
        ends = set()
        possible = set()   # the structural ends of the longest routes' leaves that can have happened
        disproved = []     # why a stated structural end of such a leaf cannot have happened
        index = _arm_index(armo)
        lacking = []       # a stated structural end whose evidence (the leaf's events) is not in the answer
        for p, f in following:
            if len(f) == best:
                pe = _leaf_ends(armo, p)
                ends |= pe
                try:
                    bad = _structural_end_problems(sc, ans, armo, p, index, arm)
                except Missing as m:   # an end without its evidence shows none of its ends either
                    lacking.append(str(m))
                    continue
                if bad:   # a leaf with an end that cannot have happened shows none of its ends
                    disproved += bad
                else:
                    possible |= pe & _STATED_STRUCTURAL_ENDS
        # a resource end needs an arm a cap stopped (G2 of the verification of lane S3)
        unstopped = _resource_end_problem(r, armo, side, ends) or next(
            (m for m in (_resource_end_below_frontier(armo, side, p) for p, f in following if len(f) == best) if m),
            None)
        gap = _computed_block(sc, ans, r, i, t, copies[0], side, L, arm) if exact is not None else None
        odd_base = _outside_alphabet(frames, side, following, best)
        if ends & _STATED_STRUCTURAL_ENDS and exact is None and possible:
            notes.append("%s: the record followed for %d bp, ended at a stated %s (%d expected otherwise)" % (
                side, best, sorted(possible)[0], expected))
        elif ends & _STATED_STRUCTURAL_ENDS and exact is None and mismatch:
            # the walk ran another strategy than the scenario's: its stop is that strategy's, named as such
            probs.append("%s: %s; the record followed for %d of %d bp to a stated %s" % (
                side, mismatch[0], best, expected, sorted(ends & _STATED_STRUCTURAL_ENDS)[0]))
        elif ends & _STATED_STRUCTURAL_ENDS and exact is None and lacking and not disproved:
            unshown.append("%s: %s (the record followed for %d of %d bp)" % (side, lacking[0], best, expected))
        elif ends & _STATED_STRUCTURAL_ENDS and exact is None:
            probs.append("a stated stop that cannot have happened: %s: %s; the record followed for %d of %d bp" % (
                side, disproved[0], best, expected))
        elif gap is not None and gap["d"] == best and _gap_stated(gap, ends, armo, following, best):
            # the record's own continuation is inadmissible there (§6.6), which the truth's minimum did not count
            searched.append("%s: the record followed for %d of %d bp, ended by %s, which the record and the seed "
                            "compute there: a truth gap" % (side, best, expected, gap["reason"]))
        elif unstopped:
            probs.append("a stated stop that cannot have happened: %s, the record followed for %d of %d bp" % (
                unstopped, best, expected))
        elif ends & _RESOURCE_ENDS or (armo.get("status") == "truncated" and trig in _RESOURCE_ENDS and
                                       not (_finite(armo.get("complete_to_bp")) and best < armo["complete_to_bp"])):
            # the cap stopped the record's route only where it reaches the depth the cap left incomplete: a path
            # that ended before complete_to_bp ended for its own stated reason (SPEC-labeled-traversal-core §6.10)
            stopped.append("%s: the record followed for %d of %d bp, stopped by %s" % (
                side, best, expected, sorted(ends & _RESOURCE_ENDS)[0] if ends & _RESOURCE_ENDS else trig))
        elif odd_base:
            # the record's next base is no base of the graph's alphabet ($ACGT, SPEC-pattern-search §8.2): no k-mer
            # holds it, so the record's route ends there, which the truth's minimum did not count
            searched.append("%s: the record followed for %d of %d bp, ended by %s before the record's base %r, "
                            "which no k-mer of the graph holds: a truth gap" % (
                                side, best, expected, sorted(ends) or "nothing stated", odd_base))
        else:
            probs.append("%s: the record followed for only %d of %d bp, the path ended by %s%s" % (
                side, best, expected, sorted(ends) or "nothing stated", how))
    if a.get("expect_branching") and not branching and not treeless:   # an arm not shown may hold the split
        bys = [((t.get("walk") or {}).get(x) or {}).get("branch_by_bp") for x in a.get("arms") or _arms_of(sc, i)]
        if stopped:
            stopped.append("no split or branch event before the stop (the record branches by %s bp)" % bys)
        else:
            probs.append("no split or branch event stated on a repeated seed (the record branches by %s bp)" % bys)
    if stopped:
        for side in a.get("arms") or _arms_of(sc, i):
            m = _arm_stop_problem(sc, ans, r, i, side)
            if m:
                probs.append("a stated stop that cannot have happened: %s" % m)
    if probs:
        return FAIL, "; ".join((probs + ["field not in the answer: %s" % x for x in unshown])[:6])
    if unshown:
        # fail on a direct arm, not exposed on the MCP's (evaluate's Missing rule); raised only where nothing
        # failed, so that a failure is never hidden behind it
        raise Missing("; ".join(unshown))
    if stopped or searched:
        return INCONCLUSIVE, "; ".join(stopped + searched + notes)
    return PASS, "; ".join(notes)


def _outside_alphabet(frames: Sequence[Tuple[Callable[[int], str], int, int]], side: str,
                      following: Sequence[Tuple[dict, str]], best: int) -> Optional[str]:
    """The record's base right after a longest route (|best| bp) along one of |frames| when it is no base the graph
    holds (not A, C, G or T: an IUPAC code of the record), else None."""
    longest = set(f for _, f in following if len(f) == best)
    for cont, room, _ in frames:
        if room <= best or cont(best) not in longest:
            continue
        nxt = cont(best + 1)
        base = (nxt[-1] if side == "right" else nxt[0]).upper() if nxt else ""
        if base and base not in "ACGT":
            return base
    return None


def _first_branch_bp(armo: dict) -> Optional[int]:
    """The depth of an arm's first split or branch event (a split's at_bp, or its segment's end in a graphlet
    arm; an event's at_bp), None without one."""
    segs = {x.get("id"): x for x in armo.get("segments") or []}
    at = []
    for sp in armo.get("splits") or []:
        if isinstance(sp.get("at_bp"), int):
            at.append(sp["at_bp"])
        elif sp.get("segment") in segs:
            sg = segs[sp["segment"]]
            at.append(sg.get("from_bp", 0) + sg.get("length_bp", 0))
    at += [e["at_bp"] for e in armo.get("branch_events") or [] if isinstance(e.get("at_bp"), int)]
    return min(at) if at else None


def check_paths_end(sc, ans, a, arm):
    """Every path of every requested arm ends with one of args.reasons_in (its path_reason, an end reason or the
    reason of a label_end event on its leaf at its end: _leaf_ends), at most args.max_bp long when given (a semantic
    stop such as below_min_labels, SPEC-labeled-traversal-core §6.4: the deployed engine states it as the end reason
    `branch` with the text below_min_labels in the label_end event, §13 I3 note 2; a graphlet as the R record's
    qualifier, `Rb`)."""
    i = a.get("seed", 0)
    r = _result(ans, i)
    if "arms" not in r:
        return FAIL, "seed %d not walked (%r)" % (i, r.get("error"))
    n = 0
    for side in _arms_of(sc, i):
        armo = r["arms"].get(side) or {}
        paths = armo.get("paths")
        if paths is None:
            raise Missing("arms.%s.paths" % side)
        if not paths:
            return FAIL, "%s: no path" % side
        for p in paths:
            ends = _leaf_ends(armo, p)
            if not ends & set(a["reasons_in"]):
                return FAIL, "%s path %s ended by %s, expected %s" % (side, p.get("id"), sorted(x for x in ends if x),
                                                                    a["reasons_in"])
            if "max_bp" in a and (p.get("length_bp") or 0) > a["max_bp"]:
                return FAIL, "%s path %s of %s bp, at most %d expected" % (side, p.get("id"), p.get("length_bp"), a["max_bp"])
            n += 1
    return PASS, "%d path(s) ended by %s" % (n, a["reasons_in"])


def check_labels_per_node(sc, ans, a, arm):
    """labels.max_labels_per_node echoed on every requested arm (labels_per_node.cap == args.cap), and a cut list
    (nodes_truncated > 0) stated by a limitation naming the knob (SPEC-labeled-traversal-core §7.0)."""
    i = a.get("seed", 0)
    r = _result(ans, i)
    if "arms" not in r:
        return FAIL, "seed %d not walked (%r)" % (i, r.get("error"))
    notes = []
    for side in _arms_of(sc, i):
        armo = r["arms"].get(side) or {}
        lpn = armo.get("labels_per_node")
        if not isinstance(lpn, dict):
            raise Missing("arms.%s.labels_per_node" % side)
        if lpn.get("cap") != a["cap"]:
            return FAIL, "%s: labels_per_node.cap %r, the request's max_labels_per_node %r" % (side, lpn.get("cap"), a["cap"])
        if (lpn.get("nodes_truncated") or 0) > 0 and not [x for x in armo.get("limitations") or []
                                                         if x.get("knob") == "labels.max_labels_per_node"]:
            return FAIL, "%s: %s node lists cut, no limitation names labels.max_labels_per_node" % (
                side, lpn.get("nodes_truncated"))
        notes.append("%s cap %s, max_seen %s, %s cut" % (side, lpn.get("cap"), lpn.get("max_seen"), lpn.get("nodes_truncated")))
    return PASS, "; ".join(notes)


def check_label_dict_contains(sc, ans, a, arm):
    """The seed's label dictionary names args.name (with args.kind when given): a derived permitted set of the
    kind asked (labels.seed_label_kind) holds the record's own label."""
    i = a.get("seed", 0)
    r = _result(ans, i)
    ld = r.get("label_dict")
    if ld is None:
        raise Missing("results[%d].label_dict" % i)
    hits = [x for x in ld if x.get("name") == a["name"] and ("kind" not in a or x.get("kind") == a["kind"])]
    if not hits:
        return FAIL, "%s (%s) not in the label dictionary of %d label(s) (%s)" % (
            a["name"], a.get("kind", "any kind"), len(ld), [(x.get("name"), x.get("kind")) for x in ld[:4]])
    others = sorted({x.get("kind") for x in ld} - ({a["kind"]} if "kind" in a else set()))
    if "kind" in a and a.get("only_kind") and others:
        return FAIL, "labels of kind %s beside %s" % (others, a["kind"])
    return PASS, "%s (%s) among %d label(s)" % (a["name"], hits[0].get("kind"), len(ld))


def check_blocked_event(sc, ans, a, arm):
    """Every requested arm of a walked seed states a blocked successor of reason args.reason (a `blocked` event,
    SPEC-labeled-traversal-core §6.6: a step into a seed node is reported as `rejoined_seed`, an edge already used
    as `edge_reuse`) at a depth of at most args.max_at_bp. A tandem-repeat seed's period closes a cycle inside the
    seed, so each arm's first step into the next copy re-enters the seed at 0 bp; the depth-0 level is processed
    whatever the budget (§6.8), so no stated stop can hide it."""
    i = a.get("seed", 0)
    r = _result(ans, i)
    if "arms" not in r:
        return FAIL, "seed %d not walked (%r)" % (i, r.get("error"))
    notes = []
    for side in a.get("arms") or _arms_of(sc, i):
        armo = (r.get("arms") or {}).get(side)
        if armo is None:
            return FAIL, "arm %s absent" % side
        segs = armo.get("segments")
        if segs is None:
            raise Missing("arms.%s.segments" % side)
        hits = [ev for sg in segs for ev in _seg_events(sg)
                if ev.get("type") == "blocked" and ev.get("reason") == a["reason"]
                and isinstance(ev.get("at_bp"), int) and ev["at_bp"] <= a.get("max_at_bp", 0)]
        if not hits:
            seen = sorted({(ev.get("reason"), ev.get("at_bp")) for sg in segs for ev in _seg_events(sg)
                           if ev.get("type") == "blocked"}, key=str)
            return FAIL, "%s: no blocked successor %s at or before %d bp (blocked events %s)" % (
                side, a["reason"], a.get("max_at_bp", 0), seen[:6])
        notes.append("%s: %s at %d bp" % (side, a["reason"], min(ev["at_bp"] for ev in hits)))
    return PASS, "; ".join(notes)


def check_walk_end(sc, ans, a, arm):
    """Near the record's end the arm stops exactly at the record's end (header label), with a structural or
    label end reason: a route that follows the record (_record_routes: any parent at a merge) reaches its end, and
    its leaf's end reasons (end_reasons and label_end events, _leaf_ends) name one of args.reasons_in. Where the
    truth's flank does not reach the record's end, the arm's leaves are read by their lengths alone. An arm a stated
    resource limit stopped before it (status truncated, every path ended by a resource reason short of the end)
    cannot show the end: inconclusive, naming the stop; so is a route search that stopped at its bound."""
    i = a.get("seed", 0)
    r = _result(ans, i)
    if "arms" not in r:
        return FAIL, "seed not walked (%r)" % r.get("error")
    rejected = _mgt_rejected(r)
    if rejected:
        return FAIL, rejected
    t = sc["truth"]
    side = a["arm"]
    armo = r["arms"][side]
    want = t["record_length"] - (t["start"] + t["length"]) if side == "right" else t["start"]
    fl = _arm_flanks(armo, side)
    lens = sorted({len(f) for _, f in fl})
    frames = [_seed_frame(t, {"start": t["start"], "flank": t["flank"], "flank_start": t["flank_start"]}, side,
                          t["length"])]
    capped = False
    labels, seq_only = _route_labels(r, t, armo)
    if frames[0][1] >= want:
        routes, _, capped = _record_routes(armo, side, frames, labels=labels)
        at_end = [(p, f) for p, f in routes if len(f) == want]
    else:
        at_end = [(p, f) for p, f in fl if len(f) == want]
    if not at_end:
        if capped:
            return INCONCLUSIVE, "the route search stopped at %d visits before a route following the record reached " \
                                 "its end (%s arm path lengths %s)" % (_MAX_ROUTE_STATES, side, lens[:8])
        if fl and armo.get("status") == "truncated" and all(
                len(f) < want and (p.get("path_reason") in _RESOURCE_ENDS or set(p.get("end_reasons") or {}) & _RESOURCE_ENDS)
                for p, f in fl):
            trig = armo.get("cap_trigger") or {}
            # a cap ends walks at or past complete_to_bp (SPEC-labeled-traversal-core §6.10), as walk_in_record
            # judges a resource end
            bad = next((m for m in (_resource_end_problem(r, armo, side, {p.get("path_reason")} | set(
                p.get("end_reasons") or {})) or _resource_end_below_frontier(armo, side, p) for p, _ in fl) if m),
                       None) or _arm_stop_problem(sc, ans, r, i, side)
            if bad:
                return FAIL, "not reached: a stated stop that cannot have happened: %s" % bad
            return INCONCLUSIVE, "not reached: the %s arm stopped by %s at %s bp (a stated stop), the record ends %d bp out" % (
                side, trig.get("reason"), trig.get("at_bp"), want)
        if want in lens:
            return FAIL, "%s arm: a path of %d bp reaches the record's end, but no route through the segments follows " \
                         "the record there" % (side, want)
        return FAIL, "%s arm path lengths %s, the record ends %d bp from the seed" % (side, lens[:8], want)
    for p, f in at_end:
        reasons = _leaf_ends(armo, p, path_reason=False)
        if a.get("reasons_in") and not reasons & set(a["reasons_in"]):
            return FAIL, "path at the record's end ended by %s, expected one of %s" % (sorted(reasons), a["reasons_in"])
    return PASS, "%s arm ends at the record's end (%d bp)%s" % (side, want, " [%s]" % seq_only if seq_only else "")


def check_traverse_seed_outcome(sc, ans, a, arm):
    r = _result(ans, a.get("seed", 0))
    walks = (r.get("outcome") or {}).get("walks")
    wants = a["walks"] if isinstance(a["walks"], list) else [a["walks"]]
    if walks not in wants:
        return FAIL, "outcome.walks %r, expected %r (error %r)" % (walks, wants, r.get("error"))
    return PASS, "outcome.walks %r" % walks


def _limitations_naming(r: dict, knob: str) -> List[dict]:
    out = [x for x in r.get("limitations") or [] if x.get("knob") == knob]
    for armo in (r.get("arms") or {}).values():
        out += [x for x in armo.get("limitations") or [] if x.get("knob") == knob]
    return out


def check_traverse_stop(sc, ans, a, arm):
    """A stated stop of the knob that was set (SPEC-labeled-traversal-core §7.0): outcome.walks among args.walks;
    with args.resource (time, memory, work: a request budget) a resource_stop of that resource whose `requested` is
    the request's knob value; with args.cap_reason, on a walked seed some requested arm truncated (or pruned) with
    that cap_trigger reason and every other one either so or complete (§6.7: an arm whose walks all ended for a
    semantic reason before the cap stays complete, for a per-arm cap and a per-seed stop alike); and a limitation
    naming args.knob (a walk_domain on the arm, or on a failed seed)."""
    i = a.get("seed", 0)
    r = _result(ans, i)
    walks = (r.get("outcome") or {}).get("walks")
    wants = a.get("walks") or ["partial"]
    if walks not in wants:
        return FAIL, "outcome.walks %r, expected %r (error %r)" % (walks, wants, r.get("error"))
    knob = a["knob"]
    bounds = ((sc["request"].get("strategy") or {}).get("bounds") or {})
    if a.get("resource"):
        rs = r.get("resource_stop")
        if not isinstance(rs, dict):
            return FAIL, "no resource_stop for a %s budget (outcome %s)" % (a["resource"], walks)
        if rs.get("resource") != a["resource"]:
            return FAIL, "resource_stop.resource %r, expected %r" % (rs.get("resource"), a["resource"])
        want = _number(bounds.get(knob.split(".", 1)[1]))
        if want is not None and _number(rs.get("requested")) != want:
            return FAIL, "resource_stop.requested %r, the request's %s is %r" % (rs.get("requested"), knob, want)
    if a.get("cap_reason") and "arms" in r:
        capped, ended = [], []
        for side in _arms_of(sc, i):
            armo = r["arms"].get(side) or {}
            trig = (armo.get("cap_trigger") or {}).get("reason")
            if armo.get("status") == "complete":
                ended.append(side)
            elif trig == a["cap_reason"] and armo.get("status") in ("truncated", "pruned"):
                capped.append(side)
            else:
                return FAIL, "%s arm status %r, cap_trigger %r, expected truncated by %s or complete" % (
                    side, armo.get("status"), trig, a["cap_reason"])
        if not capped:
            return FAIL, "no requested arm truncated by %s (%s complete)" % (a["cap_reason"], "/".join(ended))
    lims = _limitations_naming(r, knob)
    if not lims:
        return FAIL, "no limitation names %s" % knob
    return PASS, "walks %s, stopped by %s%s (%s)" % (
        walks, knob, (" (resource %s)" % a["resource"]) if a.get("resource") else "",
        ", ".join(sorted({str(x.get("kind")) for x in lims})))


def check_trace_coordinates(sc, ans, a, arm):
    """With support trace and coordinates: the seed's occurrence in the record is [start, start + length), and
    each arm's run of the record's label is adjacent to it (SPEC §7.1, record coordinates)."""
    r = _result(ans, a.get("seed", 0))
    t = sc["truth"]
    co = r.get("coordinates")
    if co is None:
        if "coordinates" not in r:
            raise Missing("coordinates")
        return FAIL, "coordinates null (%r)" % r.get("coordinates_reason")
    dict_names = [d.get("name") for d in r.get("label_dict", [])]
    if t["record"] not in dict_names:
        return FAIL, "record %s not in label_dict %r" % (t["record"], dict_names[:5])
    lid = dict_names.index(t["record"])
    s, e = t["start"], t["start"] + t["length"]
    seed_occ = [x for x in co.get("seed", []) if x.get("label") == lid]
    if not seed_occ or [s, e] not in seed_occ[0].get("occurrences", []):
        return FAIL, "seed occurrence [%d, %d) not stated (%s)" % (s, e, seed_occ[:1])
    for side, runs in (co.get("arms") or {}).items():
        for run in runs:
            if run.get("label") != lid or run.get("from_bp") != 0:
                continue
            L = run["to_bp"] - run["from_bp"]
            want = [e, e + L] if side == "right" else [s - L, s]
            if want not in run.get("occurrences", []):
                return FAIL, "%s run of %d bp not at %s (%s)" % (side, L, want, run.get("occurrences", [])[:3])
    return PASS, "seed at [%d, %d) and the arms adjacent" % (s, e)


def _seed_identity(seed: Any) -> str:
    """A seed's identity as the MCP service merges seeds (app/traversal/planner.py normalize_seeds: the sequence
    cleaned as resolve.extract_sequence cleans it, whitespace removed and upper-cased, a single FASTA record's header
    dropped; the explicit labels trimmed and de-duplicated in order, None without labels; a seed_id is no part of
    it)."""
    if not isinstance(seed, dict):
        return json.dumps(seed, sort_keys=True)
    seq = seed.get("sequence")
    if isinstance(seq, str):
        lines = [x.strip() for x in seq.splitlines()]
        seq = "".join("".join(x for x in lines if not x.startswith(">")).split()).upper()
    labels = seed.get("labels")
    if isinstance(labels, list):
        labels = list(dict.fromkeys(x.strip() if isinstance(x, str) else x for x in labels))
    return json.dumps([seq, labels], sort_keys=True)


def _distinct_seeds(seeds: Sequence[Any]) -> int:
    """The number of distinct seeds of a request by _seed_identity."""
    return len(set(_seed_identity(x) for x in seeds or []))


def _traverse_invariants(sc, ans, a, arm):
    if ans.get("http") != 200:
        return PASS, "no answer to check (http %s)" % ans.get("http")
    body = _body(ans)
    res = _get(body, "results")
    seeds = sc["request"].get("seeds", [])
    probs: List[str] = []
    unshown: List[str] = []   # the structural ends whose evidence the answer does not hold (Missing)
    if len(res) != len(seeds):
        distinct = _distinct_seeds(seeds)
        # B5b-07: a request with identical seeds states how many are distinct ("N results for M seeds (D distinct)"):
        # the service's merge onto the first gives exactly N == D (the table's traverse.duplicate_seeds explains that
        # only); a lost or a surplus result is no merge
        probs.append("%d results for %d seeds%s" % (len(res), len(seeds), " (%d distinct)" % distinct
                                                    if distinct < len(seeds) else ""))
    # the MCP walk ran the strategy the scenario declares (§5): a changed or dropped preset knob is named
    probs += _preset_mismatch(sc, ans, arm)
    # a clamp states only what the engine can lower (_clamp_problem): never a request's knob it rewrote
    probs += _clamp_problems(sc, ans, arm)
    maxbp = ((sc["request"].get("strategy") or {}).get("bounds") or {}).get("max_extension_bp")
    echoed = ((body.get("strategy") or {}).get("bounds") or {}).get("max_extension_bp")
    radius = echoed if isinstance(echoed, int) else maxbp
    for i in range(len(res)):
        r = _result(ans, i)
        if "outcome" not in r:
            probs.append("result %d without outcome" % i)
            continue
        if "arms" not in r:
            continue
        walks = (r.get("outcome") or {}).get("walks")
        for side in _arms_of(sc, i):
            armo = r["arms"].get(side)
            if armo is None:
                probs.append("result %d: arm %s absent" % (i, side))
                continue
            if armo.get("status") not in ("complete", "truncated", "pruned"):
                probs.append("result %d %s: status %r" % (i, side, armo.get("status")))
            # SPEC §7.0: complete_to_bp equals max_extension_bp exactly when the arm is complete, and walks
            # "complete" only when every requested arm is
            if isinstance(radius, int) and "complete_to_bp" in armo:
                if (armo.get("status") == "complete") != (armo.get("complete_to_bp") == radius):
                    probs.append("result %d %s: status %s with complete_to_bp %s of max_extension_bp %s" % (
                        i, side, armo.get("status"), armo.get("complete_to_bp"), radius))
            if walks == "complete" and armo.get("status") != "complete":
                probs.append("result %d: walks complete but the %s arm is %s" % (i, side, armo.get("status")))
            # the engine writes `events`, a list of event objects, on every segment of a tree or full detail
            # (traverse.cpp arm_to_json; never truncated), and graphlet_arms lists them on every segment it parses:
            # a segment without the list holds back the evidence of the decisions on it, on every arm (an MCP answer
            # with JSON segments included), whatever ends its paths claim
            bare = [s.get("id") for s in armo.get("segments") or [] if isinstance(s, dict) and
                    (not isinstance(s.get("events"), list) or not all(isinstance(x, dict) for x in s["events"]))]
            if bare:
                probs.append("result %d %s: segment(s) %s without an events list of objects (the engine writes one "
                             "on every segment of a tree or full detail)" % (i, side, bare[:6]))
            segs = {s["id"]: s for s in armo.get("segments", [])}
            for p in armo.get("paths", []):
                chain = p.get("segments", [])
                total = sum(segs[s]["length_bp"] for s in chain if s in segs)
                if total != p.get("length_bp"):
                    probs.append("result %d %s path %s: segments sum %d != length_bp %s" % (i, side, p.get("id"), total, p.get("length_bp")))
                for x, y in zip(chain, chain[1:]):
                    if y in segs and x not in (segs[y].get("parents") or []):
                        probs.append("result %d %s path %s: segment %s not a child of %s" % (i, side, p.get("id"), y, x))
                        break
                if maxbp and p.get("length_bp", 0) > maxbp:
                    probs.append("result %d %s path %s: %s bp above max_extension_bp %s" % (i, side, p.get("id"), p.get("length_bp"), maxbp))
        # a path ended by a resource reason lies on an arm that cap stopped (§7.0)
        probs += _path_resource_problems(sc, r, i)
        # a structural end lies where the request, the strategy and the arm let it happen (§6.4, §6.6)
        for side in _arms_of(sc, i):
            armo = r["arms"].get(side)
            index = _arm_index(armo) if isinstance(armo, dict) else None
            for p in (armo.get("paths") or []) if isinstance(armo, dict) else []:
                try:
                    probs += ["result %d %s %s" % (i, side, m)
                              for m in _structural_end_problems(sc, ans, armo, p, index, arm)]
                except Missing as m:
                    unshown.append("result %d %s: %s" % (i, side, m))
    probs += _traverse_reason_problems(sc, ans)
    if probs:
        return FAIL, "; ".join((probs + ["field not in the answer: %s" % x for x in unshown])[:6])
    if unshown:
        # an end without its evidence: fail on a direct arm, not exposed on the MCP's (evaluate's Missing rule),
        # raised only where nothing else failed, so that a failure is never hidden behind it
        raise Missing(unshown[0])
    return PASS, "%d result(s) consistent" % len(res)


def check_invariants(sc, ans, a, arm, client=None):
    route = sc.get("route", "")
    if route == "POST /pattern":
        return _pattern_invariants(sc, ans, a, arm)
    if route == "POST /resolve":
        return _resolve_invariants(sc, ans, a, arm, client)
    if route == "POST /traverse":
        return _traverse_invariants(sc, ans, a, arm)
    return PASS, "no invariants for %s" % route


# ---------------------------------------------------------------------------------------------- registry
# The checks that hold on any answer of the route, whatever the scenario asks: the HTTP status, the route's
# invariants, the counts' consistency and the time budget. Every other check is a PROPERTY check, the claim the
# scenario was built to test; an execution whose property checks all could not test anything (not exposed on the
# MCP arm) is inconclusive even when its generic checks pass (scenario_outcome). A scenario with no property check
# (a refusal: its status, code and text are the claim) is judged by its generic checks.
GENERIC_CHECKS = frozenset({"status", "invariants", "counts_consistent", "budget_respected"})


def is_property(name: str) -> bool:
    """Whether check |name| tests the scenario's own claim (README, "Outcomes")."""
    return name not in GENERIC_CHECKS


CHECKS: Dict[str, Callable] = {
    "status": check_status,
    "answer_fields": check_answer_fields,
    "entry_fields": check_entry_fields,
    "deterministic": check_deterministic,
    "content_encoding": check_content_encoding,
    "capabilities_doc": check_capabilities_doc,
    "pattern_slot_error": check_pattern_slot_error,
    "count_relation": check_count_relation,
    "withheld": check_withheld,
    "cut": check_cut,
    "stop": check_stop,
    "retrieval_complete": check_retrieval_complete,
    "notes": check_notes,
    "low_complexity_note": check_low_complexity_note,
    "stops_stated": check_stops_stated,
    "results_levels": check_results_levels,
    "absent": check_absent,
    "contains_source_context": check_contains_source_context,
    "source_absent": check_source_absent,
    "source_label": check_source_label,
    "source_placement": check_source_placement,
    "path_spells": check_path_spells,
    "path_record_verified": check_path_record_verified,
    "path_absent": check_path_absent,
    "predicate_selected": check_predicate_selected,
    "predicate_either_mirror": check_predicate_either_mirror,
    "motif": check_motif,
    "clamped": check_clamped,
    "budget_respected": check_budget_respected,
    "counts_consistent": check_counts_consistent,
    "invariants": check_invariants,
    "resolve_runs": check_resolve_runs,
    "label_not_full": check_label_not_full,
    "resolve_absent": check_resolve_absent,
    "resolve_graph_runs": check_resolve_graph_runs,
    "resolve_discover_contains": check_resolve_discover_contains,
    "walk_in_record": check_walk_in_record,
    "walk_end": check_walk_end,
    "blocked_event": check_blocked_event,
    "traverse_seed_outcome": check_traverse_seed_outcome,
    "traverse_stop": check_traverse_stop,
    "labels_per_node": check_labels_per_node,
    "paths_end": check_paths_end,
    "label_dict_contains": check_label_dict_contains,
    "trace_coordinates": check_trace_coordinates,
}


# ---------------------------------------------------------------------------------------------- relations
# The manifest's `relations` (README.md, "Relations"): comparisons of two scenarios' answers on the same arm. The
# comparator (lane S2b's compare.py) evaluates them; this is the reference implementation of their semantics.


def _drop(o: Any, ignore: Sequence[str], prefix: str = "") -> Any:
    """|o| without every `timing` key and without the dotted paths in |ignore|."""
    if isinstance(o, dict):
        out = {}
        for k, v in o.items():
            p = prefix + k
            if k == "timing" or p in ignore:
                continue
            out[k] = _drop(v, ignore, p + ".")
        return out
    if isinstance(o, list):
        return [_drop(x, ignore, prefix) for x in o]
    return o


def _json_leaf_equal(a: Any, b: Any) -> bool:
    """Two JSON leaves equal as JSON values (suite_common.json_equal's leaf rule): a boolean only the same boolean
    (true is not 1, false not 0; the owner's decision on SD-011), numbers by value (1 equals 1.0)."""
    if isinstance(a, bool) or isinstance(b, bool):
        return isinstance(a, bool) and isinstance(b, bool) and a == b
    return a == b


def _first_difference(a: Any, b: Any, path: str = "") -> str:
    """The first difference of two JSON values as text, or "": objects key by key, lists item by item, leaves as JSON
    values (_json_leaf_equal); compare.py's first_difference is the same rule (tests/test_compare_report.py)."""
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b)):
            if k not in a or k not in b:
                return "%s%s only in %s" % (path, k, "a" if k in a else "b")
            d = _first_difference(a[k], b[k], path + k + ".")
            if d:
                return d
        return ""
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            return "%s: %d vs %d items" % (path.rstrip(".") or "<root>", len(a), len(b))
        for i, (x, y) in enumerate(zip(a, b)):
            d = _first_difference(x, y, "%s%d." % (path, i))
            if d:
                return d
        return ""
    return "" if _json_leaf_equal(a, b) else "%s: %s vs %s" % (path.rstrip(".") or "<root>", json.dumps(a)[:80],
                                                             json.dumps(b)[:80])


def _relation_same_answer(rel, a_sc, A, b_sc, B):
    args = rel.get("args") or {}
    if A.get("http") != B.get("http"):
        return FAIL, "http %s vs %s" % (A.get("http"), B.get("http"))
    if A.get("http") != 200:
        ea, eb = A.get("error") or {}, B.get("error") or {}
        if (ea.get("code"), ea.get("text")) != (eb.get("code"), eb.get("text")):
            return FAIL, "errors differ: %r vs %r" % (ea, eb)
        return PASS, "the same refusal (http %s)" % A.get("http")
    ba, bb = _body(A), _body(B)
    if args.get("paths"):
        for pa, pb in args["paths"]:
            va = ba if pa == "" else _get(ba, pa)
            vb = bb if pb == "" else _get(bb, pb)
            d = _first_difference(_drop(va, ()), _drop(vb, ()))
            if d:
                return FAIL, "a.%s differs from b.%s: %s" % (pa or "<root>", pb or "<root>", d)
        return PASS, "%d field(s) equal" % len(args["paths"])
    ignore = args.get("ignore") or []
    d = _first_difference(_drop(ba, ignore), _drop(bb, ignore))
    if d:
        return FAIL, "answers differ: %s" % d
    return PASS, "equal answers (timing%s aside)" % ((", " + ", ".join(ignore)) if ignore else "")


def _kmer_set(runs: Any, offset: int, lo: int, hi: int) -> set:
    out = set()
    for r in runs or []:
        for x in range(max(r[0] + offset, lo), min(r[1] + offset, hi)):
            out.add(x)
    return out


def _expects_stop(sc: dict) -> bool:
    """Whether a /resolve scenario's property is a stated stop: its resolve_runs asks prefix_if_stopped (a budget
    built to stop inside the query)."""
    return any(isinstance(e, dict) and e.get("check") == "resolve_runs" and (e.get("args") or {}).get("prefix_if_stopped")
               for e in sc.get("expect") or [])


def _relation_prefix_equals(rel, a_sc, A, b_sc, B):
    outcome, detail = _prefix_comparison(rel, A, B)
    if outcome == PASS and _expects_stop(a_sc) and not (_body(A) or {}).get("stop"):
        # a's property is its stop: two complete profiles that agree say nothing of the prefix a stop answers
        return INCONCLUSIVE, "no stop happened on a: the prefix property was not tested (%s)" % detail
    return outcome, detail


def _prefix_comparison(rel, A, B):
    if A.get("http") != 200 or B.get("http") != 200:
        return INCONCLUSIVE, "http %s / %s: no two profiles to compare" % (A.get("http"), B.get("http"))
    args = rel.get("args") or {}
    ba, bb = _body(A), _body(B)
    fa, fb = int(args.get("a_from_kmer", 0)), int(args.get("b_from_kmer", 0))
    na, nb = _get(ba, "num_kmers"), _get(bb, "num_kmers")
    lo, hi = max(fa, fb), min(fa + na, fb + nb)
    if hi <= lo:
        return INCONCLUSIVE, "no common k-mer (a [%d, %d), b [%d, %d))" % (fa, fa + na, fb, fb + nb)
    ga, gb = _kmer_set(_get(ba, "graph_runs"), fa, lo, hi), _kmer_set(_get(bb, "graph_runs"), fb, lo, hi)
    if ga != gb:
        return FAIL, "graph runs differ on k-mers [%d, %d): %d vs %d in the graph (first difference %s)" % (
            lo, hi, len(ga), len(gb), min(ga ^ gb))
    for name in args.get("labels") or []:
        la, lb = _resolve_label(ba, name), _resolve_label(bb, name)
        sa, sb = _kmer_set(la.get("runs"), fa, lo, hi), _kmer_set(lb.get("runs"), fb, lo, hi)
        if sa != sb:
            return FAIL, "label %s differs on k-mers [%d, %d): %d vs %d supported (first difference %s)" % (
                name, lo, hi, len(sa), len(sb), min(sa ^ sb))
    stops = ["%s stopped at %s" % (n, x["stop"].get("resolved_kmers")) for n, x in (("a", ba), ("b", bb)) if x.get("stop")]
    return PASS, "agree on k-mers [%d, %d)%s" % (lo, hi, (" (" + ", ".join(stops) + ")") if stops else "")


def _relation_zero_equals_absent(rel, a_sc, A, b_sc, B):
    if A.get("http") != 200 or B.get("http") != 200:
        return INCONCLUSIVE, "http %s / %s" % (A.get("http"), B.get("http"))
    args = rel.get("args") or {}
    i = int(args.get("pattern", 0))
    e = _entry(A, i)
    if "error" in e:
        return INCONCLUSIVE, "pattern %d refused (%s)" % (i, e["error"].get("code"))
    value, relation_ = _count(_get(e, "counts." + args.get("count", "contexts")))
    if relation_ != "exact":
        return INCONCLUSIVE, "count %s %s: not an absence claim either way" % (relation_, value)
    _, pat = _request_pattern(a_sc, i)
    n = len(pat) - K + 1
    full = _get(_body(B), "graph_runs") == [[0, n]]
    if (value == 0) == (not full):
        return PASS, "count exact %d and %s" % (value, "every k-mer in the graph" if full else "a k-mer absent")
    return FAIL, "count exact %d but graph_runs %s (of %d k-mers)" % (value, _get(_body(B), "graph_runs")[:4], n)


def _result_keys(e: dict, key: str, level: Optional[int]) -> List[tuple]:
    res = e.get("results")
    if res is None:
        raise Missing("results")
    out = []
    for r in res:
        if level is not None and r.get("mismatches", 0) != level:
            continue
        st = r.get("strand", r.get("orientation"))
        out.append((r.get("sequence"), st) if key == "path" else (r.get("kmer"), r.get("offset"), st))
    return out


def _level_count(e: dict, level: Optional[int]) -> Optional[Tuple[Optional[int], str]]:
    c = (e.get("counts") or {}).get("contexts") or (e.get("counts") or {}).get("paths")
    if not isinstance(c, dict):
        return None
    if level is not None:
        c = (c.get("by_mismatches") or {}).get(str(level))
        if not isinstance(c, dict):
            return None
    return c.get("value"), c.get("relation")


def _relation_subset(rel, a_sc, A, b_sc, B):
    if A.get("http") != 200 or B.get("http") != 200:
        return INCONCLUSIVE, "http %s / %s" % (A.get("http"), B.get("http"))
    args = rel.get("args") or {}
    ea, eb = _entry(A, int(args.get("a_pattern", 0))), _entry(B, int(args.get("b_pattern", 0)))
    if "error" in ea or "error" in eb:
        return INCONCLUSIVE, "a refused slot (%s / %s)" % ((ea.get("error") or {}).get("code"), (eb.get("error") or {}).get("code"))
    key = args.get("key", "context")
    ka = _result_keys(ea, key, args.get("a_level"))
    kb = set(_result_keys(eb, key, args.get("b_level")))
    missing = [x for x in ka if x not in kb]
    if missing:
        if eb.get("retrieval_complete") is True:
            return FAIL, "%d of a's %d results not in b's complete list (first %s)" % (len(missing), len(ka), missing[0])
        return INCONCLUSIVE, "%d of a's %d results not in b's list, which is not complete (withheld %r, cut %r)" % (
            len(missing), len(ka), eb.get("withheld"), eb.get("cut"))
    ca, cb = _level_count(ea, args.get("a_level")), _level_count(eb, args.get("b_level"))
    if ca and cb and ca[1] == "exact" and cb[1] == "exact" and ca[0] > cb[0]:
        return FAIL, "a's exact count %d above b's exact count %d" % (ca[0], cb[0])
    return PASS, "a's %d results within b's %d" % (len(ka), len(kb))


RELATIONS: Dict[str, Callable] = {
    "same_answer": _relation_same_answer,
    "prefix_equals": _relation_prefix_equals,
    "zero_equals_absent": _relation_zero_equals_absent,
    "subset": _relation_subset,
}


def evaluate_relation(rel: dict, a_scenario: dict, a_answer: dict, b_scenario: dict, b_answer: dict, arm: str) -> dict:
    """A relation of the manifest on one arm: {"relation", "type", "a", "b", "outcome": pass | fail | inconclusive
    | not_exposed, "detail"}, from the two normalised answers (the first repeat of each, by convention)."""
    fn = RELATIONS.get(rel.get("type"))
    out = {"relation": rel.get("id"), "type": rel.get("type"), "a": rel.get("a"), "b": rel.get("b")}
    if fn is None:
        out.update(outcome=FAIL, detail="unknown relation type %r" % rel.get("type"))
        return out
    try:
        outcome, detail = fn(rel, a_scenario, a_answer, b_scenario, b_answer)
    except Missing as m:
        outcome, detail = (NOT_EXPOSED if arm == "mcp" else FAIL), "field not in the answer: %s" % m
    except Exception as ex:   # never take the comparison down
        outcome, detail = FAIL, "relation raised %s: %s" % (type(ex).__name__, ex)
    out.update(outcome=outcome, detail=detail)
    return out


# The arguments each check reads (required, optional); `any_status` (evaluate) is accepted by every check. The
# builder validates every expectation against this table, so that a misspelt argument cannot weaken a check.
CHECK_ARGS: Dict[str, Tuple[Tuple[str, ...], Tuple[str, ...]]] = {
    "status": ((), ("http", "code", "text", "text_prefix", "text_contains")),
    "answer_fields": ((), ("equals", "present", "absent", "when_http")),
    "entry_fields": ((), ("pattern", "equals", "present", "absent", "tolerance")),
    "deterministic": ((), ()),
    "content_encoding": (("value",), ()),
    "capabilities_doc": ((), ("required", "equals", "max_bytes")),
    "pattern_slot_error": (("pattern", "code"), ("text", "text_prefix", "text_contains", "fields")),
    "count_relation": ((), ("pattern", "count", "level", "relation", "value", "min", "max")),
    "withheld": (("reason",), ("pattern",)),
    "cut": (("reason",), ("pattern",)),
    "stop": ((), ("pattern", "null", "phase", "reason", "mismatch_level_in")),
    "retrieval_complete": (("value",), ("pattern",)),
    "notes": ((), ("pattern", "has", "lacks")),
    "low_complexity_note": (("expected",), ("pattern",)),
    "stops_stated": (("reasons",), ()),
    "results_levels": (("max_level",), ("pattern", "min_returned")),
    "absent": ((), ("pattern", "count")),
    "contains_source_context": ((), ("pattern", "truth", "strand", "level", "require", "substitutions")),
    "source_absent": ((), ("pattern", "strand")),
    "source_label": ((), ("pattern", "strand")),
    "source_placement": ((), ("pattern", "strand")),
    "path_spells": ((), ("pattern", "strand")),
    "path_record_verified": ((), ("pattern", "strand")),
    "path_absent": ((), ("pattern", "sequence")),
    "predicate_selected": (("selected",), ("pattern", "strand", "selection_strand_in", "label")),
    "predicate_either_mirror": ((), ("pattern", "label", "min")),
    "motif": ((), ("pattern", "selected", "decided_by_in", "label_present")),
    "clamped": (("field", "requested", "effective"), ("where",)),
    "budget_respected": ((), ("rule", "cap_ms")),
    "counts_consistent": ((), ("pattern",)),
    "invariants": ((), ()),
    "resolve_runs": (("label",), ("full", "prefix_if_stopped", "unless_stopped", "covers", "run_end_in",
                                  "run_start_in", "runs", "kmers", "max_kmers", "min_kmers")),
    "label_not_full": (("label",), ()),
    "resolve_absent": ((), ()),
    "resolve_graph_runs": (("runs",), ()),
    "resolve_discover_contains": (("label",), ()),
    "walk_in_record": ((), ("seed", "arms", "label_scope", "any_copy", "expect_branching", "prefix_only", "if_walked")),
    "walk_end": (("arm",), ("seed", "reasons_in")),
    "blocked_event": (("reason",), ("seed", "arms", "max_at_bp")),
    "traverse_seed_outcome": (("walks",), ("seed",)),
    "traverse_stop": (("knob",), ("seed", "walks", "resource", "cap_reason")),
    "labels_per_node": (("cap",), ("seed",)),
    "paths_end": (("reasons_in",), ("seed", "max_bp")),
    "label_dict_contains": (("name",), ("seed", "kind", "only_kind")),
    "trace_coordinates": ((), ("seed",)),
}


# the checks of one seed of a /traverse answer (their `seed` argument)
_SEED_CHECKS = frozenset(n for n, (req, opt) in CHECK_ARGS.items() if "seed" in req + opt)


def validate_expectation(exp: dict) -> List[str]:
    """What is wrong with one expectation: an unknown check, a missing or unknown argument, or an argument set
    that asserts nothing (a count_relation without a value, a minimum of at least 1, a maximum or a relation that
    excludes something)."""
    name, args = exp.get("check"), exp.get("args") or {}
    if name not in CHECKS or name not in CHECK_ARGS:
        return ["unknown check %r" % name]
    req, opt = CHECK_ARGS[name]
    probs = ["%s: argument %s missing" % (name, x) for x in req if x not in args]
    probs += ["%s: unknown argument %s" % (name, x) for x in args if x not in req and x not in opt and x != "any_status"]
    if name == "count_relation":
        rels = args.get("relation")
        rels = set(rels if isinstance(rels, list) else [rels]) if rels is not None else None
        if args.get("min") is not None and args["min"] < 1:
            probs.append("count_relation: min %r asserts nothing" % args["min"])
        if "value" not in args and "max" not in args and not args.get("min") and \
                (rels is None or rels >= {"exact", "at_least", "bounds", "unknown"}):
            probs.append("count_relation: no value, minimum, maximum or relation to assert")
    return probs


# The host's statements that an attempt's duration bound stopped a seed (SPEC-labeled-traversal-core §6.8): a stop the
# host states and enforces itself, after which nothing of that seed is delivered, so no traversal property of that
# seed can be tested. Each is matched in full, in the wording the deployed host serves and in the current source's
# (which adds the in_ram load term), never as a substring: the MCP tool's own hint ("the attempt reached the duration
# bound the server enforces for it; a retry needs a new attempt_id") is no host statement.
# AttemptAtBound (src/cli/traverse_attempts.cpp Attempt::check_delivery): the bound reached while the response was
# built or written; the host answers the request 503 `deadline` (server.cpp), which the MCP service relays as the
# seed result's error. Accepted on the MCP arm only: a direct arm receives the 503 itself, so a 200 seed result
# stating it is not the host's answer.
_AT_BOUND_RE = re.compile(
    r"the attempt reached the duration bound the server enforces for it \((\d+) ms: the "
    r"(?:in_ram load \(if any\) plus the )?seeds' time budgets plus the server's allowance, at most hard_cap_ms, the "
    r"content timeout less one second; see usage\.bound\) while its response was built or written: nothing of it "
    r"is delivered")
# attempt_deadline (src/cli/traverse.cpp external_cause and not_started_seed_to_json, src/graph/traversal/walker.cpp
# seed_external_stop): the attempt's walk-until (its bound less the larger of half the allowance and the delivery
# reserve) passed before the seed began or in its seed phase; the host answers 200 with the seed failed
# (resource_stop attempt_deadline). Any arm.
_WALK_UNTIL_CAUSE = (r"the attempt reached the time at which the server stops walking it \(the duration bound it "
                     r"enforces, (?:in_ram load \+ seeds' time budgets \+ allowance|the seeds' time budgets plus its "
                     r"allowance) \(at most hard_cap_ms\), less the larger of half the allowance and the delivery "
                     r"reserve, kept for the delivery; see usage\.bound\.walk_until_ms\)")
_NOT_STARTED_RE = re.compile(r"not started: " + _WALK_UNTIL_CAUSE + r" before this seed began, (\d+) ms after the "
                             r"request was received; no budget of the request ran out: a new attempt can traverse it")
_SEED_PHASE_STOP_RE = re.compile(
    r"the attempt reached the time at which the server stops walking it \(its duration bound less the larger of "
    r"half the allowance and the delivery reserve; see usage\.bound\.walk_until_ms\) while the seed was "
    r"(?:read to derive its permitted set|validated against its labels), (\d+) ms after the request was received: "
    r"no traversal was made")


def _seed_error_text(r: Any) -> Optional[str]:
    err = r.get("error") if isinstance(r, dict) else None
    if isinstance(err, dict):
        err = err.get("text")
    return err if isinstance(err, str) else None


def _attempt_stops(scenario: dict, answer: dict, arm: str, client: Optional[dict] = None
                   ) -> Dict[int, Tuple[str, Optional[str]]]:
    """{seed index: (the host's statement, why it cannot have happened or None)} for the seed results of a /traverse
    200 that state that the attempt's bound stopped them (above). The stated bound must have been reached where the
    answer states the attempt's elapsed time: AttemptAtBound's bound against usage.elapsed_ms and the client's wall
    time (`client`: the runner's, which spans the host's work; the MCP runner passes the service's backend time); a
    walk-until stop at or after usage.bound.walk_until_ms where the answer states it, with its resource_stop (where
    stated) naming attempt_deadline and the seed failed."""
    out: Dict[int, Tuple[str, Optional[str]]] = {}
    if not str(scenario.get("route", "")).endswith("/traverse") or answer.get("http") != 200:
        return out
    body = answer.get("body")
    if not isinstance(body, dict):
        return out
    usage = body.get("usage") if isinstance(body.get("usage"), dict) else {}
    bound_block = usage.get("bound") if isinstance(usage.get("bound"), dict) else {}
    for i, r in enumerate(body.get("results") or []):
        text = _seed_error_text(r)
        if not text or "arms" in r:
            continue
        m = _AT_BOUND_RE.fullmatch(text)
        if m:
            if arm != "mcp":
                continue
            bound = float(m.group(1))
            seen = [("usage.elapsed_ms", _number(usage.get("elapsed_ms"))),
                    ("the client's wall time", _number((client or {}).get("wall_ms")))]
            early = [(what, v) for what, v in seen if v is not None and v < bound - 1.0]
            out[i] = (text, "seed %d: the attempt's bound of %.0f ms stated as reached after %.0f ms (%s)" % (
                i, bound, early[0][1], early[0][0]) if early else None)
            continue
        m = _NOT_STARTED_RE.fullmatch(text) or _SEED_PHASE_STOP_RE.fullmatch(text)
        if not m:
            continue
        used = float(m.group(1))
        rs = r.get("resource_stop") if isinstance(r.get("resource_stop"), dict) else None
        walks = (r.get("outcome") or {}).get("walks") if isinstance(r.get("outcome"), dict) else None
        until = _number(bound_block.get("walk_until_ms"))
        problem = None
        if rs is not None and rs.get("resource") != "attempt_deadline":
            problem = "seed %d: stopped by the attempt's walk-until, but its resource_stop names %r" % (
                i, rs.get("resource"))
        elif walks is not None and walks != "failed":
            problem = "seed %d: stopped by the attempt's walk-until before any traversal, but outcome.walks %r" % (
                i, walks)
        elif until is not None and used < until - 1.0:
            problem = "seed %d: stopped by the attempt's walk-until %.0f ms after the request was received, before " \
                      "usage.bound.walk_until_ms %.0f" % (i, used, until)
        out[i] = (text, problem)
    return out


def _with_stopped_seeds_failed(answer: dict, stops: Dict[int, Tuple[str, Optional[str]]]) -> dict:
    """|answer| with each stopped seed's result in the shape of a failed seed (an outcome, no arms), so that the checks
    of the whole answer judge the other seeds and pass over these."""
    body = dict(answer["body"])
    res = list(body.get("results") or [])
    for i in stops:
        if isinstance(res[i], dict) and "outcome" not in res[i]:
            res[i] = dict(res[i], outcome={"walks": "failed"})
    body["results"] = res
    return dict(answer, body=body)


# the checks that read the runner's `client` ({"wall_ms"}) beside the answer: budget_respected (an allowed 503
# deadline, _budget_503), the /resolve stop's time (_resolve_stop_problem) and status (the front ingress's 504 at
# the cut, ingress_cut_allowed); deterministic reads it too
_CLIENT_CHECKS = {"budget_respected", "invariants", "resolve_runs", "status"}


def evaluate(scenario: dict, answer: dict, arm: str, previous: Optional[Sequence[Any]] = None,
             client: Optional[dict] = None) -> List[dict]:
    """Every expectation of |scenario| against |answer| (DESIGN §8.2), each result {check, outcome, detail,
    property} (is_property). A check whose field is missing answers not_exposed on the MCP arm and fail on a direct
    arm; an unknown check name or a check that raises is a fail stating why (never a silent pass). `deterministic` is
    left out of the first repeat (nothing to compare). `previous`: the earlier repeats' digests or {digest,
    time_limited}; `client` (optional, from a runner): {"wall_ms"} of this execution, which budget_respected reads
    on an allowed 503 `deadline` (_budget_503), status on the front ingress's 504 (ingress_cut_allowed), and the
    checks that judge a /resolve stop's time against it (_resolve_stop_problem: invariants, resolve_runs,
    deterministic)."""
    out = []
    http = answer.get("http")
    allowed = [200]
    for exp in scenario.get("expect", []):
        if exp.get("check") == "status":
            w = exp.get("args", {}).get("http", 200)
            allowed = w if isinstance(w, list) else [w]
    if deadline_allowed(scenario, answer):
        allowed = list(allowed) + [503]   # a 503 deadline: the other checks cannot test their property
    if ingress_cut_allowed(scenario, answer, client, arm):
        allowed = list(allowed) + [504]   # the front ingress's cut at the budget: nothing of the host to test
    # the seeds the attempt's duration bound stopped (_attempt_stops): their own checks are inconclusive (or fail
    # when the stated bound cannot have been reached), the other seeds keep theirs
    stops = _attempt_stops(scenario, answer, arm, client)
    stopped_view = _with_stopped_seeds_failed(answer, stops) if stops else answer
    stop_problems = [p for _, p in stops.values() if p]
    for exp in scenario.get("expect", []):
        name = exp.get("check")
        args = exp.get("args", {})
        if name == "deterministic" and not previous:
            continue   # the first repeat has nothing to compare: no result, so that it neither passes nor holds up
        fn = CHECKS.get(name)
        if fn is None:
            out.append({"check": name, "outcome": FAIL, "detail": "unknown check %r" % name, "property": True})
            continue
        if name not in ("status", "deterministic") and http != 200 and not args.get("any_status"):
            if http in allowed and name == "budget_respected" and http == 503:
                outcome, detail = _budget_503(scenario, answer, client, args)
                out.append({"check": name, "outcome": outcome, "detail": detail, "property": is_property(name)})
            elif http in allowed:
                out.append({"check": name, "outcome": INCONCLUSIVE,
                            "detail": "not testable: http %s, an answer the scenario allows (%s)" % (
                                http, ((answer.get("error") or {}).get("text") or "")[:200]),
                            "property": is_property(name)})
            else:
                out.append({"check": name, "outcome": FAIL,
                            "detail": "no answer to check: http %s (%s)" % (
                                http, ((answer.get("error") or {}).get("text") or "")[:200]),
                            "property": is_property(name)})
            continue
        judged = answer
        if stops and name not in ("status", "deterministic"):
            seed = args.get("seed", 0)
            if name in _SEED_CHECKS and seed in stops:
                text, problem = stops[seed]
                if problem:
                    out.append({"check": name, "outcome": FAIL, "property": is_property(name),
                                "detail": "a stated attempt-bound stop that cannot have happened: %s" % problem})
                else:
                    out.append({"check": name, "outcome": INCONCLUSIVE, "property": is_property(name),
                                "detail": "not testable: seed %d stopped by its attempt's duration bound (stated by "
                                          "the host: %s)" % (seed, text[:160])})
                continue
            if name not in _SEED_CHECKS:
                judged = stopped_view
        try:
            if name == "deterministic":
                outcome, detail = fn(scenario, answer, args, arm, previous=previous, client=client)
            elif name in _CLIENT_CHECKS:
                outcome, detail = fn(scenario, judged, args, arm, client=client)
            else:
                outcome, detail = fn(scenario, judged, args, arm)
            if name == "invariants" and stop_problems:
                outcome = FAIL
                detail = "a stated attempt-bound stop that cannot have happened: %s; %s" % (
                    "; ".join(stop_problems), detail)
        except Missing as m:
            outcome = NOT_EXPOSED if arm == "mcp" else FAIL
            detail = "field not in the answer: %s" % m
        except Undecidable as u:
            outcome, detail = INCONCLUSIVE, "not testable: %s" % u
        except Exception as ex:   # a check must never take the run down
            outcome, detail = FAIL, "check raised %s: %s" % (type(ex).__name__, ex)
        out.append({"check": name, "outcome": outcome, "detail": detail, "property": is_property(name)})
    return out


def scenario_outcome(results: List[dict]) -> str:
    """The execution's outcome from its checks, by precedence fail > finding > inconclusive > pass (the rule lane
    S2's runner, comparator and report and the async session's MCP runner apply): `fail` when a check failed (or
    reported an outcome outside the vocabulary); else `finding` when a check is a finding (the host broke a stated
    bound without failing the request); else `inconclusive` when a check could not test its property (an allowed
    503, a list a legitimate statement hides, a legitimate resource stop before the property), when there is no
    check or every check is not exposed, or when the scenario's PROPERTY checks (is_property) all could not test
    anything although its generic checks (status, invariants, counts_consistent, budget_respected) passed; else
    `pass`."""
    outs = [r.get("outcome") for r in results]
    if any(o == FAIL or o not in OUTCOMES for o in outs):
        return FAIL
    if FINDING in outs:
        return FINDING
    if INCONCLUSIVE in outs:
        return INCONCLUSIVE
    if not any(o != NOT_EXPOSED for o in outs):
        return INCONCLUSIVE
    props = [r for r in results if is_property(r.get("check"))]
    if props and not any(r.get("outcome") == PASS for r in props):
        return INCONCLUSIVE
    return PASS
