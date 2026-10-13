"""The suite's checks (DESIGN §8.2): pure functions of a scenario and a normalised answer, stdlib only, shared by
every arm so that a disagreement between arms is never a difference of checking code.

    evaluate(scenario, answer, arm, previous=None)
        -> [{"check": name, "outcome": "pass" | "fail" | "not_exposed" | "finding", "detail": str}]

`answer` is the NORMALISED answer: {"http": int, "body": obj | None, "error": {"code", "text"} | None} and, where
the arm has them, "headers" (lower-case names). `normalise_direct` builds it from a direct HTTP answer; the MCP
arm's normaliser maps its tool results into the same shape for the fields both expose. A check that needs a field
the MCP arm does not expose answers "not_exposed" on that arm (and "fail" on a direct arm, where the field is part
of the contract). A `finding` is a stated bound or SPEC promise the host broke without failing the request (the
time budget, budget_respected): the report lists findings in their own section. `previous` (optional) lists the
`answer_digest`s of the earlier repeats of the same scenario and arm, for the `deterministic` check.

Every check named in a manifest's `expect` is a function registered in CHECKS; new checks are added here only.
Coordinates: truth["start"] is 0-based on the record's forward strand; host placements (`nt_coords`) are 1-based
inclusive; resolve runs and traverse intervals 0-based half-open (the SPECs' conventions).
"""

from __future__ import annotations

import gzip
import hashlib
import json
import math
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
_ORIENT_RANK = {"+": 0, "forward": 0, "-": 1, "reverse": 1, "=": 2, "palindromic": 2}

PASS, FAIL, NOT_EXPOSED = "pass", "fail", "not_exposed"
# the host broke a stated bound or a SPEC promise without failing the request: reported in its own section of the
# report, never counted as a pass (budget_respected)
FINDING = "finding"


class Missing(Exception):
    """A field the check needs is absent from the answer."""


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


def expected_substitutions(kind: str, pattern: str, instance: str, strand: str, table: int = 1) -> List[dict]:
    """The substitutions the host must state for |instance| against |pattern| (SPEC-pattern-search.md §21.8):
    positions (residues) in the coordinates of P; for '-' the instance is read in P's orientation."""
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
def check_status(sc, ans, a, arm):
    http = ans.get("http")
    want = a.get("http", 200)
    wants = want if isinstance(want, list) else [want]
    if http not in wants:
        err = ans.get("error") or {}
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


def check_answer_fields(sc, ans, a, arm):
    """Fields of the answer's body by dotted path (any_status: also on an error body, e.g. a 409's state)."""
    body = _body(ans)
    for path, want in (a.get("equals") or {}).items():
        got = _get(body, path)
        if got != want:
            return FAIL, "%s = %r, expected %r" % (path, got, want)
    for path in a.get("present", []):
        _get(body, path)
    for path in a.get("absent", []):
        try:
            _get(body, path)
            return FAIL, "%s present, expected absent" % path
        except Missing:
            pass
    return PASS, "fields as expected"


def check_entry_fields(sc, ans, a, arm):
    e = _entry(ans, a.get("pattern", 0))
    for path, want in (a.get("equals") or {}).items():
        got = _get(e, path)
        if isinstance(want, float) and isinstance(got, (int, float)):
            if abs(got - want) > a.get("tolerance", 0.051):
                return FAIL, "%s = %r, expected %r" % (path, got, want)
        elif got != want:
            return FAIL, "%s = %r, expected %r" % (path, got, want)
    for path in a.get("present", []):
        _get(e, path)
    for path in a.get("absent", []):
        try:
            _get(e, path)
            return FAIL, "%s present, expected absent" % path
        except Missing:
            pass
    return PASS, "entry fields as expected"


def check_deterministic(sc, ans, a, arm, previous=None):
    d = answer_digest(ans, sc.get("route"))
    if not previous:
        return PASS, "first repeat (digest %s)" % d[:16]
    diff = [p for p in previous if p != d]
    if diff:
        return FAIL, "digest %s differs from %d earlier repeat(s) (%s)" % (d[:16], len(diff), diff[0][:16])
    return PASS, "digest %s equal to %d earlier repeat(s)" % (d[:16], len(previous))


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


def check_count_relation(sc, ans, a, arm):
    e = _entry(ans, a.get("pattern", 0))
    if "error" in e:
        return FAIL, "entry refused (%s)" % e["error"].get("code")
    name = a.get("count", "contexts")
    c = _get(e, "counts." + name)
    if "level" in a:
        c = _get(c, "by_mismatches.%d" % a["level"])
    value, rel = _count(c)
    if "relation" in a:
        rels = a["relation"] if isinstance(a["relation"], list) else [a["relation"]]
        if rel not in rels:
            return FAIL, "%s relation %s (value %s), expected %s" % (name, rel, value, rels)
    lower = c.get("lower", value) if rel == "bounds" else value
    if "value" in a and (rel != "exact" or value != a["value"]):
        return FAIL, "%s = %s %s, expected exact %s" % (name, rel, value, a["value"])
    if "min" in a:
        if rel == "unknown":
            return FAIL, "%s unknown, expected at least %d" % (name, a["min"])
        if rel == "exact" and value < a["min"]:
            return FAIL, "%s exact %d below %d" % (name, value, a["min"])
        if rel == "bounds" and c.get("upper", 0) < a["min"]:
            return FAIL, "%s bounds upper %s below %d" % (name, c.get("upper"), a["min"])
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
        if a.get("optional"):
            return PASS, "no stop (allowed)"
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
        return FAIL, "%d of %d source contexts (strand %s, level %s) returned; missing offsets %s (returned %d, %s)" % (
            len(found), len(wanted), strand, level, missing[:8], len(e.get("results", [])),
            "complete" if e.get("retrieval_complete") else "incomplete")
    if require == "any" and not found:
        return FAIL, "none of %d source contexts returned (strand %s, level %s)" % (len(wanted), strand, level)
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
    if e.get("retrieval_complete") is not True:
        return FAIL, "inconclusive: the list is not complete (withheld %r, cut %r)" % (e.get("withheld"), e.get("cut"))
    return PASS, "none of %d source contexts among %d results (complete)" % (len(wanted), len(e.get("results", [])))


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
        return FAIL, "no source context returned"
    for r in found:
        if r.get("labels") is None:
            raise Missing("results[].labels")
        if _label_of(r, t["taxid"]):
            return PASS, "taxid %s among the labels of the source context (offset %d)" % (t["taxid"], r["offset"])
    return FAIL, "taxid %s on none of %d source contexts (labels_status %s)" % (
        t["taxid"], len(found), [r.get("labels_status") for r in found])


def check_source_placement(sc, ans, a, arm):
    """The record is placed: the source context's label (the record's taxid) lists an occurrence of the record at
    the window's 1-based coordinates and strand (SPEC §14.3)."""
    i = a.get("pattern", 0)
    e = _entry(ans, i)
    t = sc["truth"]
    strand = a.get("strand", t.get("strand", "+"))
    wanted, found = _find_source(e, t, strand, None, sc["request"].get("scope", "any_offset"))
    if not found:
        return FAIL, "no source context returned"
    coords = "%d-%d" % (t["start"] + 1, t["start"] + t["length"])
    for r in found:
        lab = _label_of(r, t["taxid"])
        if lab is None:
            continue
        occ = lab.get("occurrence_list")
        if occ is None:
            raise Missing("occurrence_list")
        for o in occ:
            if o.get("record") == t["record"] and o.get("nt_coords") == coords and o.get("strand") == strand:
                if "seq_id" in o and o["seq_id"] != t["seq_id"]:
                    return FAIL, "record %s placed with seq_id %s, the index's is %s" % (t["record"], o["seq_id"], t["seq_id"])
                if o.get("nt_length") not in (None, t["record_length"]):
                    return FAIL, "nt_length %s, the record's length is %s" % (o.get("nt_length"), t["record_length"])
                return PASS, "%s %s %s placed (offset %d)" % (t["record"], coords, strand, r["offset"])
    return FAIL, "no occurrence %s %s %s under taxid %s on %d source contexts" % (
        t["record"], coords, strand, t["taxid"], len(found))


def _paths_with_sequence(e: dict, seq: str) -> List[dict]:
    res = e.get("results")
    if res is None:
        raise Missing("results")
    return [r for r in res if r.get("sequence") == seq]


def check_path_spells(sc, ans, a, arm):
    e = _entry(ans, a.get("pattern", 0))
    t = sc["truth"]
    strand = a.get("strand", t.get("strand", "+"))
    hits = [r for r in _paths_with_sequence(e, t["window"]) if r.get("strand", r.get("orientation")) in (
        [strand] if strand != "=" else ["=", "palindromic"])]
    if not hits:
        return FAIL, "no path spells the record's bases (strand %s; %d paths returned, %s)" % (
            strand, len(e.get("results", [])), "complete" if e.get("retrieval_complete") else "incomplete")
    for r in hits:
        if r.get("instance") not in (None, r["sequence"]) or r.get("offset", 0) != 0:
            return FAIL, "path instance/offset inconsistent"
        n = t["length"] - K + 1
        if "nodes" in r and len(r["nodes"]) != n:
            return FAIL, "path of %d nodes, expected %d" % (len(r["nodes"]), n)
    return PASS, "%d path(s) spell the record's bases" % len(hits)


def check_path_record_verified(sc, ans, a, arm):
    e = _entry(ans, a.get("pattern", 0))
    t = sc["truth"]
    strand = a.get("strand", t.get("strand", "+"))
    coords = "%d-%d" % (t["start"] + 1, t["start"] + t["length"])
    hits = _paths_with_sequence(e, t["window"])
    if not hits:
        return FAIL, "no path spells the record's bases (%d paths, %s)" % (
            len(e.get("results", [])), "complete" if e.get("retrieval_complete") else "incomplete")
    for r in hits:
        if r.get("labels") is None:
            continue
        lab = _label_of(r, t["taxid"])
        if lab is None:
            continue
        if lab.get("support") != "record_verified":
            return FAIL, "taxid %s carries the path with support %r" % (t["taxid"], lab.get("support"))
        for o in lab.get("occurrence_list") or []:
            if o.get("record") == t["record"] and o.get("nt_coords") == coords and o.get("strand") == strand:
                return PASS, "%s %s %s record_verified" % (t["record"], coords, strand)
        return FAIL, "taxid %s record_verified, but no occurrence %s %s %s (%s)" % (
            t["taxid"], t["record"], coords, strand, (lab.get("occurrence_list") or [])[:3])
    if all(r.get("labels") is None for r in hits):
        raise Missing("path labels")
    return FAIL, "taxid %s on none of the %d path(s) spelling the record's bases" % (t["taxid"], len(hits))


def check_path_absent(sc, ans, a, arm):
    e = _entry(ans, a.get("pattern", 0))
    seq = a.get("sequence") or sc["truth"].get("spliced_dna")
    hits = _paths_with_sequence(e, seq)
    if hits:
        return FAIL, "%d path(s) spell the sequence expected absent" % len(hits)
    if e.get("retrieval_complete") is not True:
        return FAIL, "inconclusive: the list is not complete (withheld %r, cut %r)" % (e.get("withheld"), e.get("cut"))
    return PASS, "no path spells it (complete list of %d)" % len(e.get("results", []))


def check_predicate_selected(sc, ans, a, arm):
    """The source contexts are selected (in the results) or not (absent from a complete filtered list)."""
    e = _entry(ans, a.get("pattern", 0))
    if "error" in e:
        return FAIL, "entry refused (%s)" % e["error"].get("code")
    t = sc["truth"]
    strand = a.get("strand", t.get("strand", "+"))
    wanted, found = _find_source(e, t, strand, None, sc["request"].get("scope", "any_offset"))
    if a["selected"]:
        if len(found) != len(wanted):
            return FAIL, "%d of %d source contexts selected" % (len(found), len(wanted))
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
    if e.get("retrieval_complete") is not True:
        return FAIL, "inconclusive: the filtered list is not complete"
    return PASS, "source contexts not selected (complete filtered list)"


def check_predicate_either_mirror(sc, ans, a, arm):
    """Under predicate_strands "either": every returned '-' context whose reverse complement is a record k-mer
    of the truth is selected by the label with selection strand reverse_complement or both (SPEC §19.5)."""
    e = _entry(ans, a.get("pattern", 0))
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
            if r["selection_strands"][labs.index(col)] not in ("reverse_complement", "both"):
                return FAIL, "mirror context selected with strand %r" % r["selection_strands"][labs.index(col)]
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


def check_budget_respected(sc, ans, a, arm):
    """The server's elapsed time against the request's effective time budget, by each route's SPEC rule, with the
    allowance args.allowance_ms taken from the SPEC and the capabilities snapshot (never fitted to an answer).
    /pattern and /resolve: the top-level timing.elapsed_ms within limits.time_budget_ms (allowance 0: past the
    budget their SPECs answer 503, so a 200 was written within it). /traverse (SPEC-labeled-traversal-core §6.8):
    each walked seed's deadline is its budget on the seed's clock, which never stops the seed phase (validation or
    derivation, depth 0); so the walk must end within the allowance (the walk-stop time from a deadline to the
    walk's end) plus the walk's longest uninterruptible piece after max(budget, seed_phase_ms), and a seed phase
    that outlasted the budget leaves one level (complete_to_bp <= 1). A breach is a FINDING (the request still
    answered); a seed phase longer than the budget is a timing fact stated in the detail, not a finding."""
    allowance = _number(a.get("allowance_ms")) or 0.0
    body = _body(ans)
    route = sc.get("route", "")
    req = sc.get("request") if isinstance(sc.get("request"), dict) else {}
    if route == "POST /traverse":
        bounds = (req.get("strategy") or {}).get("bounds") or {}
        budget = _number(((body.get("strategy") or {}).get("bounds") or {}).get("time_budget_ms"))
        if budget is None:
            budget = _number(bounds.get("time_budget_ms"))
        if budget is None:
            raise Missing("strategy.bounds.time_budget_ms")
        res = _get(body, "results")
        over, facts, seen = [], [], 0
        for i, r in enumerate(res):
            t = r.get("timing")
            if not isinstance(t, dict) or _number(t.get("elapsed_ms")) is None:
                if "arms" in r:
                    raise Missing("results[%d].timing.elapsed_ms" % i)
                continue   # a failed seed without its timing: nothing to compare
            if "arms" not in r:
                continue   # no walk: the seed phase alone, which the budget does not govern
            seen += 1
            el = _number(t["elapsed_ms"])
            sp = _number(t.get("seed_phase_ms"))
            if sp is None:
                raise Missing("results[%d].timing.seed_phase_ms" % i)
            dl = t.get("deadline") or {}
            lp = dl.get("longest_piece") or {}
            lp_ms = _number(lp.get("ms")) or 0.0
            # a longest piece longer than the walk phase lay in the seed phase: it grants the walk nothing
            piece = lp_ms if lp_ms <= el - sp else 0.0
            late = el - max(budget, sp)
            evidence = "elapsed %.0f ms, seed phase %.0f ms, longest piece %s %.0f ms (%s rows), stopped_by %s" % (
                el, sp, lp.get("kind"), lp_ms, lp.get("rows"), dl.get("stopped_by"))
            if late > allowance + piece:
                over.append("seed %d: the walk ended %.0f ms after max(budget %.0f, seed phase %.0f), above the "
                            "allowance %.0f + its longest piece %.0f (%s)" % (i, late, budget, sp, allowance, piece, evidence))
            if sp >= budget:
                facts.append("seed %d: the seed phase (%.0f ms) outlasted the budget, which does not govern it" % (i, sp))
                for side in _arms_of(sc, i):
                    armo = (r.get("arms") or {}).get(side) or {}
                    if (armo.get("cap_trigger") or {}).get("reason") == "time_budget" and \
                            (armo.get("complete_to_bp") or 0) > 1:
                        over.append("seed %d %s: the budget passed in the seed phase, yet the walk went on to "
                                    "complete_to_bp %s (one level expected)" % (i, side, armo.get("complete_to_bp")))
        if over:
            return FINDING, "; ".join(over[:4])
        if not seen:
            return INCONCLUSIVE, "no walked seed: nothing the budget governs"
        return PASS, "%d walked seed(s) ended within the allowance %.0f ms of max(budget %.0f, seed phase)%s" % (
            seen, allowance, budget, ("; " + "; ".join(facts)) if facts else "")
    budget = _number((body.get("limits") or {}).get("time_budget_ms"))
    if budget is None:
        named = req.get("time_budget_ms") if route == "POST /pattern" else (req.get("bounds") or {}).get("time_budget_ms")
        budget = _number(named)
    if budget is None:
        raise Missing("limits.time_budget_ms")
    el = _number(_get(body, "timing.elapsed_ms"))
    if el is None:
        raise Missing("timing.elapsed_ms")
    if el > budget + allowance:
        stops = []
        for i, e in enumerate(body.get("patterns") or []):
            if isinstance(e, dict) and e.get("stop"):
                stops.append("p%d stop %s/%s" % (i, e["stop"].get("phase"), e["stop"].get("reason")))
        if body.get("stop"):
            stops.append("stop %s/%s" % (body["stop"].get("phase"), body["stop"].get("reason")))
        return FINDING, "elapsed %.0f ms above the time budget %.0f ms + allowance %.0f ms (%s)" % (
            el, budget, allowance, ", ".join(stops) or "no stop stated")
    return PASS, "elapsed %.0f ms within %.0f + %.0f ms" % (el, budget, allowance)


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
                subs = expected_substitutions(kind, pat, inst, strand, table) if strand != "=" else None
                if subs is not None:
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
    if e.get("retrieval_complete"):
        name = _primary_count_name(e, req)
        c = counts.get(name)
        if not c or c.get("relation") != "exact" or c.get("value") != len(res):
            probs.append("p%d: retrieval_complete with %s %r and %d results" % (i, name, c and (c.get("relation"), c.get("value")), len(res)))
    stop = e.get("stop")
    if stop and stop.get("reason") in _BUDGET_STOP_REASONS and "low_complexity_pattern" in (e.get("notes") or []):
        probs.append("p%d: low_complexity_pattern beside a budget stop" % i)
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


def check_counts_consistent(sc, ans, a, arm):
    body = _body(ans)
    pats = _get(body, "patterns")
    limits = body.get("limits") or {}
    which = [a["pattern"]] if "pattern" in a else range(len(pats))
    probs: List[str] = []
    for i in which:
        probs += _counts_consistent_problems(sc, i, pats[i], limits)
    # one budget per request: a budget stop is sticky for the later patterns (SPEC §7.6)
    first = None
    for i, e in enumerate(pats):
        s = e.get("stop")
        if first is not None and "counts" in e and "no_stop_codon" not in (e.get("notes") or []):
            if not s or s.get("reason") != first.get("reason"):
                probs.append("p%d: after a %s stop of an earlier pattern, stop %r" % (i, first.get("reason"), s))
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
    limits = body.get("limits") or {}
    for i, e in enumerate(pats):
        req_id = sc["request"]["patterns"][i].get("id") if isinstance(sc["request"]["patterns"][i], dict) else None
        if e.get("id") != req_id:
            probs.append("p%d: id %r, request %r" % (i, e.get("id"), req_id))
        probs += _pattern_entry_problems(sc, i, e, limits, None)
    return (FAIL, "; ".join(probs[:6])) if probs else (PASS, "%d entries consistent" % len(pats))


# ---------------------------------------------------------------------------------------------- resolve checks
def _resolve_label(body: dict, name: str) -> dict:
    labs = [l for l in _get(body, "labels") if l.get("label") == name]
    if not labs:
        raise Missing("label %s" % name)
    return labs[0]


def check_resolve_runs(sc, ans, a, arm):
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
    if a.get("unless_stopped") and body.get("stop"):
        return PASS, "stopped after %s k-mers: the prefix's answer (invariants check it)" % body["stop"].get("resolved_kmers")
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


def check_resolve_discover_contains(sc, ans, a, arm):
    body = _body(ans)
    names = [l.get("label") for l in _get(body, "labels")]
    trunc = body.get("labels_truncated")
    if a["label"] in names:
        return PASS, "%s discovered (%d labels%s)" % (a["label"], len(names), ", truncated" if trunc else "")
    if trunc:
        return PASS, "%s not listed, the list is truncated (kept %s of %s; index order, not best-supported)" % (
            a["label"], trunc.get("kept"), trunc.get("total"))
    return FAIL, "%s not discovered in a complete list of %d labels" % (a["label"], len(names))


def _resolve_invariants(sc, ans, a, arm):
    if ans.get("http") != 200:
        return PASS, "no answer to check (http %s)" % ans.get("http")
    body = _body(ans)
    req = sc["request"]
    probs = []
    stop = body.get("stop")
    seq = req.get("sequence", "") if isinstance(req, dict) else ""
    qk = max(0, len(seq) - K + 1)
    num = body.get("num_kmers")
    if stop:
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


def graphlet_arms(mgt: str) -> Dict[str, dict]:
    """The arms of an MGT v1 body (SPEC-labeled-traversal-core.md §7.5) in the shape of a `detail: full` arm, as
    far as the checks read it: segments (id, parents, from_bp, length_bp, sequence in natural orientation),
    splits (non-empty when a segment ends in a split), branch_events (the V records), and paths (one per leaf,
    the chain through parents[0], length_bp = from_bp + length_bp, end_reasons from the runs anchored at the leaf
    with to_bp = the path length, path_reason from T)."""
    arms: Dict[str, dict] = {}
    cur = None
    seg = None
    for line in mgt.split("\n"):
        if not line:
            continue
        f = line.split(" ")
        rec = f[0]
        if rec == "A":
            side = "left" if f[1] == "l" else "right"
            cur = {"status": {"c": "complete", "t": "truncated", "p": "pruned"}.get(f[2], f[2]),
                   "segments": [], "splits": [], "branch_events": [], "paths": [], "_runs": [], "_leaves": {}}
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
            seg = {"id": len(cur["segments"]), "parents": parents, "from_bp": int(f[2]), "length_bp": int(f[3]),
                   "sequence": natural}
            if f[10] == "*":
                del seg["sequence"]
            cur["segments"].append(seg)
            if f[8] in ("0", "1"):
                cur["splits"].append({"segment": seg["id"], "ambiguous": f[8] == "1"})
        elif rec == "V":
            cur["branch_events"].append({"at_bp": int(f[1]), "segment": int(f[2])})
        elif rec == "T" and seg is not None:
            cur["_leaves"][seg["id"]] = None if f[1] == "*" else _END_CODES.get(f[1][0], f[1])
        elif rec == "R":
            cur["_runs"].append({"segment": int(f[1]), "label": int(f[2]), "from_bp": int(f[3]), "to_bp": int(f[4]),
                                 "end": f[5]})
    for side, a in arms.items():
        segs = a["segments"]
        for leaf, reason in sorted(a["_leaves"].items()):
            chain = [leaf]
            while segs[chain[-1]]["parents"]:
                chain.append(segs[chain[-1]]["parents"][0])
            chain.reverse()
            length = segs[leaf]["from_bp"] + segs[leaf]["length_bp"]
            ends: Dict[str, int] = {}
            for r in a["_runs"]:
                if r["segment"] == leaf and r["to_bp"] == length and r["end"] != "m":
                    # the label's end reason as the full result counts it: the code's enum (Dh -> dead_end)
                    name = _END_CODES.get(r["end"][0], r["end"])
                    ends[name] = ends.get(name, 0) + 1
            a["paths"].append({"id": len(a["paths"]), "segments": chain, "length_bp": length, "end_reasons": ends,
                               "path_reason": reason})
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
        arms = {}
        for side, a in parsed.items():
            merged = dict((r.get("arms") or {}).get(side) or {})
            merged.update(a)
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


def _arms_of(sc: dict, i: int) -> List[str]:
    d = (sc["request"].get("strategy") or {}).get("direction", "both")
    return ["left", "right"] if d == "both" else [d]


def check_walk_in_record(sc, ans, a, arm):
    """The walk follows the record: per requested arm, some path's flank is the record's own continuation from
    the seed (for its whole length; with any_copy, the continuation of any copy of a repeated seed), and the root
    segment, which every path shares, agrees with it. With the record's header as the only label (label_scope
    "record") no flank may run past the record's end. expect_branching: some arm states a split or a branch
    event."""
    i = a.get("seed", 0)
    r = _result(ans, i)
    if "arms" not in r:
        return FAIL, "seed %d not walked (outcome %r, error %r)" % (i, r.get("outcome"), r.get("error"))
    t = sc["truth"] if "truths" not in sc else sc["truths"][i]
    L = t["length"]
    copies = [{"start": t["start"], "flank": t["flank"], "flank_start": t["flank_start"]}]
    if a.get("any_copy"):
        copies += t.get("copies", [])

    def continuation(c: dict, side: str, n: int) -> str:
        s0 = c["start"] - c["flank_start"]
        return c["flank"][s0 + L:s0 + L + n] if side == "right" else c["flank"][max(0, s0 - n):s0]

    def room(c: dict, side: str) -> int:
        s0 = c["start"] - c["flank_start"]
        return len(c["flank"]) - (s0 + L) if side == "right" else s0

    rec_end_dist = t["record_length"] - (t["start"] + L)
    probs, notes = [], []
    branching = False
    for side in a.get("arms") or _arms_of(sc, i):
        armo = r["arms"].get(side)
        if armo is None:
            probs.append("arm %s absent" % side)
            continue
        if armo.get("splits") or armo.get("branch_events"):
            branching = True
        fl = _arm_flanks(armo, side)
        if not fl:
            notes.append("%s: no path" % side)
            continue
        follows, best, beyond = 0, 0, 0
        for p, f in fl:
            if a.get("label_scope") == "record" and len(f) > (rec_end_dist if side == "right" else t["start"]):
                probs.append("%s path %s runs %d bp, past the record's end" % (side, p.get("id"), len(f)))
                continue
            ok = False
            for c in copies:
                if room(c, side) < len(f):
                    beyond += 1
                    continue
                if f == continuation(c, side, len(f)):
                    ok = True
                    break
            if ok:
                follows += 1
                best = max(best, len(f))
        if follows == 0:
            lens = [len(f) for _, f in fl]
            probs.append("%s: none of %d paths follows the record (lengths %s%s)" % (
                side, len(fl), lens[:6], "; %d beyond the truth flank" % beyond if beyond else ""))
            continue
        for s in [s for s in armo.get("segments", []) if not s.get("parents")]:
            rs = s.get("sequence", "")
            if not any(room(c, side) >= len(rs) and rs == continuation(c, side, len(rs)) for c in copies):
                probs.append("%s: the root segment disagrees with the record" % side)
        notes.append("%s: %d of %d paths follow the record (longest %d bp)" % (side, follows, len(fl), best))
    if a.get("expect_branching") and not branching:
        probs.append("no split or branch event stated on a repeated seed")
    if probs:
        return FAIL, "; ".join(probs[:6])
    return PASS, "; ".join(notes)


def check_walk_end(sc, ans, a, arm):
    """Near the record's end the arm stops exactly at the record's end (header label), with a structural or
    label end reason. An arm a resource stopped before it (a stated stop: status truncated, every path ended by
    a resource reason short of the end) cannot show the end: the check passes saying so, the walk so far being
    walk_in_record's to check."""
    i = a.get("seed", 0)
    r = _result(ans, i)
    if "arms" not in r:
        return FAIL, "seed not walked (%r)" % r.get("error")
    t = sc["truth"]
    side = a["arm"]
    armo = r["arms"][side]
    want = t["record_length"] - (t["start"] + t["length"]) if side == "right" else t["start"]
    fl = _arm_flanks(armo, side)
    lens = sorted({len(f) for _, f in fl})
    if want not in lens:
        resource = {"max_steps", "time_budget", "resource_limit", "max_live_paths", "max_paths", "max_output_bp"}
        if fl and armo.get("status") == "truncated" and all(
                len(f) < want and (p.get("path_reason") in resource or set(p.get("end_reasons") or {}) & resource)
                for p, f in fl):
            trig = armo.get("cap_trigger") or {}
            return PASS, "not reached: the %s arm stopped by %s at %s bp (a stated stop), the record ends %d bp out" % (
                side, trig.get("reason"), trig.get("at_bp"), want)
        return FAIL, "%s arm path lengths %s, the record ends %d bp from the seed" % (side, lens[:8], want)
    for p, f in fl:
        if len(f) == want:
            reasons = set((p.get("end_reasons") or {}).keys())
            if a.get("reasons_in") and not reasons & set(a["reasons_in"]):
                return FAIL, "path at the record's end ended by %s, expected one of %s" % (sorted(reasons), a["reasons_in"])
    return PASS, "%s arm ends at the record's end (%d bp)" % (side, want)


def check_traverse_seed_outcome(sc, ans, a, arm):
    r = _result(ans, a.get("seed", 0))
    walks = (r.get("outcome") or {}).get("walks")
    wants = a["walks"] if isinstance(a["walks"], list) else [a["walks"]]
    if walks not in wants:
        return FAIL, "outcome.walks %r, expected %r (error %r)" % (walks, wants, r.get("error"))
    return PASS, "outcome.walks %r" % walks


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


def _traverse_invariants(sc, ans, a, arm):
    if ans.get("http") != 200:
        return PASS, "no answer to check (http %s)" % ans.get("http")
    body = _body(ans)
    res = _get(body, "results")
    seeds = sc["request"].get("seeds", [])
    probs = []
    if len(res) != len(seeds):
        probs.append("%d results for %d seeds" % (len(res), len(seeds)))
    maxbp = ((sc["request"].get("strategy") or {}).get("bounds") or {}).get("max_extension_bp")
    for i in range(len(res)):
        r = _result(ans, i)
        if "outcome" not in r:
            probs.append("result %d without outcome" % i)
            continue
        if "arms" not in r:
            continue
        for side in _arms_of(sc, i):
            armo = r["arms"].get(side)
            if armo is None:
                probs.append("result %d: arm %s absent" % (i, side))
                continue
            if armo.get("status") not in ("complete", "truncated", "pruned"):
                probs.append("result %d %s: status %r" % (i, side, armo.get("status")))
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
    return (FAIL, "; ".join(probs[:6])) if probs else (PASS, "%d result(s) consistent" % len(res))


def check_invariants(sc, ans, a, arm):
    route = sc.get("route", "")
    if route == "POST /pattern":
        return _pattern_invariants(sc, ans, a, arm)
    if route == "POST /resolve":
        return _resolve_invariants(sc, ans, a, arm)
    if route == "POST /traverse":
        return _traverse_invariants(sc, ans, a, arm)
    return PASS, "no invariants for %s" % route


# ---------------------------------------------------------------------------------------------- registry
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
    "traverse_seed_outcome": check_traverse_seed_outcome,
    "trace_coordinates": check_trace_coordinates,
}


# ---------------------------------------------------------------------------------------------- relations
# The manifest's `relations` (README.md, "Relations"): comparisons of two scenarios' answers on the same arm. The
# comparator (lane S2b's compare.py) evaluates them; this is the reference implementation of their semantics.
INCONCLUSIVE = "inconclusive"


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


def _first_difference(a: Any, b: Any, path: str = "") -> str:
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
    return "" if a == b else "%s: %s vs %s" % (path.rstrip(".") or "<root>", json.dumps(a)[:80], json.dumps(b)[:80])


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


def _relation_prefix_equals(rel, a_sc, A, b_sc, B):
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


def evaluate(scenario: dict, answer: dict, arm: str, previous: Optional[Sequence[str]] = None) -> List[dict]:
    """Every expectation of |scenario| against |answer| (DESIGN §8.2). A check whose field is missing answers
    not_exposed on the MCP arm and fail on a direct arm; an unknown check name or a check that raises is a fail
    stating why (never a silent pass)."""
    out = []
    http = answer.get("http")
    allowed = [200]
    for exp in scenario.get("expect", []):
        if exp.get("check") == "status":
            w = exp.get("args", {}).get("http", 200)
            allowed = w if isinstance(w, list) else [w]
    for exp in scenario.get("expect", []):
        name = exp.get("check")
        args = exp.get("args", {})
        fn = CHECKS.get(name)
        if fn is None:
            out.append({"check": name, "outcome": FAIL, "detail": "unknown check %r" % name})
            continue
        if name not in ("status", "deterministic") and http != 200 and not args.get("any_status"):
            if http in allowed:
                out.append({"check": name, "outcome": PASS,
                            "detail": "not applicable: http %s, an answer the scenario allows" % http})
            else:
                out.append({"check": name, "outcome": FAIL,
                            "detail": "no answer to check: http %s (%s)" % (
                                http, ((answer.get("error") or {}).get("text") or "")[:200])})
            continue
        try:
            if name == "deterministic":
                outcome, detail = fn(scenario, answer, args, arm, previous=previous)
            else:
                outcome, detail = fn(scenario, answer, args, arm)
        except Missing as m:
            outcome = NOT_EXPOSED if arm == "mcp" else FAIL
            detail = "field not in the answer: %s" % m
        except Exception as ex:   # a check must never take the run down
            outcome, detail = FAIL, "check raised %s: %s" % (type(ex).__name__, ex)
        out.append({"check": name, "outcome": outcome, "detail": detail})
    return out


def scenario_outcome(results: List[dict]) -> str:
    """fail when a check failed; finding when none failed and a check is a finding (the host broke a stated bound
    without failing the request: never counted as a pass); pass otherwise (not_exposed checks aside)."""
    if any(r["outcome"] == FAIL for r in results):
        return "fail"
    if any(r["outcome"] == FINDING for r in results):
        return FINDING
    return "pass"
