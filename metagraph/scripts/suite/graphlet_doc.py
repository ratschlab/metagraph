"""The documents of area G as the suite reads them (DESIGN-v8 §2 "bodies first", FORMAT-v8 §3 `doc`): the two digests
by which the body and envelope pre-checks compare what each way processed, computed by the suite's own reading of the
MGT v1 frame (SPEC-labeled-traversal-core.md §7.5.1, §7.5.2), never by the library. Stdlib only, Python 3.8+.

    frame(text)             -> {"head", "seed_index", "j", "lines", "body_start", "z_start"}   (ValueError naming why)
    body_digest(text)       -> sha256 hex of the body: H with its seed_index masked, J and Z left out
    envelope(text)          -> the J line's object, or None for a server body (no J)
    envelope_digest(text)   -> sha256 hex of the envelope without its per-request values (timing, usage, results[*]
                               .timing, results[*].annotation.direct_reads/keys_mapped, results[*].seed.seed_id and
                               .seed_id_mismatch), or None (no J)

A server body (`detail: graphlet`, results[i].graphlet) is `H S X* L* O Q? K* (A B* V* (G P* E* T? C?)* R*)* Z`; the
standalone file of the same seed (the library's save/standalone_text, the MCP's mgt export) has the J line right after
H, and its Z counts one line more. So the body digest leaves out J and Z, and both forms of one retrieval have one
digest. It masks H's seed_index (the seed's position in the request, SPEC §7.5.2: a batch position, not the
retrieval), so the same seed sent alone (W1/W2, index 0) and in a service batch has one digest too.

The library's parser.body_digest (parser.py:931-945) hashes the same bytes with seed_index NOT masked:
body_digest(text, mask_seed_index=False) equals it on every document (tests/test_graphlet_local.py); that is the one
intended difference. Unlike the library's, this reader refuses a document whose frame is broken (no final LF, no H of
MGT v1 first, a J anywhere but right after H, no Z last counting the lines): a cut body never gets a digest.
"""

from __future__ import annotations

import hashlib
import json
import re
from typing import Any, Optional

# SPEC §7.5.2: H mgt 1 <k> <regime> <alphabet> <mode> <support> <reconverge> <cap> <continuation_bp> <seed_index>
# walk <index_ns|*> <index_fp|*> <index_meta_fp>: 16 fields, seed_index the 12th (index 11)
H_FIELDS = 16
H_SEED_INDEX = 11
SEED_INDEX_MASK = "*"   # no valid H carries `*` there (an integer, §7.5.1), so a masked H is never a real one
_MGT_INT = re.compile(r"(0|[1-9][0-9]*)\Z")
# the per-request values of a J envelope (DESIGN-v8 §2: W1/W2's J equals W4's after dropping them): the response's
# timing and its attempt usage, and the seed summary's timing (SPEC §7.5.5: "timing (per seed and per response)").
# strategy.output.timing is a strategy setting, not a measurement, and stays.
ENVELOPE_DROP = ("timing", "usage")
SUMMARY_DROP = ("timing",)
# FORMAT v2 E1: the seed summary's annotation counters of lookahead work, which differ between two engine runs of one
# request (SPEC-labeled-traversal-core.md §8.3: direct_reads and keys_mapped "count lookahead work"; INVENTORY §4.1).
# The strategy's own `annotation` block (batch_kmers) is a setting and stays, as do access_path and rows_requested.
ANNOTATION_DROP = ("direct_reads", "keys_mapped")
# VERIFY-L VL16 (the coordinator's decisions): the seed summary's seed.seed_id is request metadata (the caller's own
# name for the seed; the service may send its own, the direct G requests none), and so is seed.seed_id_mismatch (true
# when a supplied seed_id differs from the recomputed one, SPEC-core §5 step 3); validated_seed_id and the rest stay
SEED_DROP = ("seed_id", "seed_id_mismatch")


def _text(text: Any) -> str:
    """|text| as a str: bytes are decoded as UTF-8 strictly (an MGT document is UTF-8, §7.5.1)."""
    if isinstance(text, (bytes, bytearray)):
        try:
            return bytes(text).decode("utf-8")
        except UnicodeDecodeError as e:
            raise ValueError("not UTF-8: %s" % e) from None
    if not isinstance(text, str):
        raise ValueError("an MGT document is text, not %s" % type(text).__name__)
    return text


def frame(text: Any) -> dict:
    """The frame of an MGT v1 document (a server body or a standalone file), checked as SPEC §7.5.1/§7.5.2 state it:
    LF-terminated records (LF is the only separator: names may hold other line breaks raw), the first an H record of
    MGT v1 with 16 fields and an integer seed_index, an optional J record right after it (a JSON object: the library's
    envelope line, never the server's), the last a Z record counting the document's lines. Nothing else is read.
    -> {"head": H's fields, "seed_index": int, "j": the J line's JSON text or None, "lines": the line count,
    "body_start": where the records after H (and J) begin, "z_start": where the Z line begins}. A ValueError names the
    rule a document breaks."""
    text = _text(text)
    if not text.endswith("\n"):
        raise ValueError("truncated document: it does not end with a line feed")
    first = text.find("\n")
    head = text[:first].split(" ")
    if len(head) != H_FIELDS or head[:3] != ["H", "mgt", "1"]:
        raise ValueError("line 1: not the H record of MGT v1 (H mgt 1 and 16 fields), %r" % text[:first][:80])
    if not _MGT_INT.match(head[H_SEED_INDEX]):
        raise ValueError("line 1: H seed_index %r is not an integer without sign or leading zero"
                         % head[H_SEED_INDEX])
    body_start = first + 1
    j = None
    if text.startswith("J ", body_start):
        end = text.find("\n", body_start)
        j = text[body_start + 2:end]
        body_start = end + 1
    lines = text.count("\n")
    z_start = text.rfind("\n", 0, len(text) - 1) + 1
    z = text[z_start:-1].split(" ")
    if z_start < body_start or len(z) != 2 or z[0] != "Z" or not _MGT_INT.match(z[1]) or int(z[1]) != lines:
        raise ValueError("line %d: not the Z record counting the document's %d lines, %r"
                         % (lines, lines, text[z_start:-1][:80]))
    # a J line is the library's and comes right after H only (§7.5.2): one anywhere else is no document of either form
    if text.find("\nJ ", body_start - 1, z_start) >= 0:
        raise ValueError("a J record after line %d: the envelope line follows H only" % (2 if j is not None else 1))
    return {"head": head, "seed_index": int(head[H_SEED_INDEX]), "j": j, "lines": lines, "body_start": body_start,
            "z_start": z_start}


def body_digest(text: Any, mask_seed_index: bool = True) -> str:
    """sha256 (hex) of the document's body as UTF-8: the H line (its seed_index field written `*` when
    |mask_seed_index|) with its LF, then every record after H and J up to, not including, the Z record. A server body
    and the standalone file saved from it have one digest; so have two retrievals of one seed at different batch
    positions. mask_seed_index=False gives the library's parser.body_digest (parser.py:931-945)."""
    text = _text(text)
    f = frame(text)
    head = list(f["head"])
    if mask_seed_index:
        head[H_SEED_INDEX] = SEED_INDEX_MASK
    h = hashlib.sha256((" ".join(head) + "\n").encode("utf-8"))
    h.update(text[f["body_start"]:f["z_start"]].encode("utf-8"))
    return h.hexdigest()


def envelope(text: Any) -> Optional[dict]:
    """The J line's object (the response envelope with results[] reduced to this seed's summary, SPEC §7.5.2), None
    for a document without J (a server body). A J line that is not a JSON object is a ValueError."""
    j = frame(text)["j"]
    if j is None:
        return None
    try:
        obj = json.loads(j)
    except ValueError as e:
        raise ValueError("line 2: the J record is not JSON: %s" % e) from None
    if not isinstance(obj, dict):
        raise ValueError("line 2: the J record is not a JSON object")
    return obj


def _summary_view(r: Any) -> Any:
    if not isinstance(r, dict):
        return r
    out = {k: v for k, v in r.items() if k not in SUMMARY_DROP}
    if isinstance(out.get("annotation"), dict):
        out["annotation"] = {k: v for k, v in out["annotation"].items() if k not in ANNOTATION_DROP}
    if isinstance(out.get("seed"), dict):
        out["seed"] = {k: v for k, v in out["seed"].items() if k not in SEED_DROP}
    return out


def envelope_view(env: dict) -> dict:
    """|env| without its per-request values (FORMAT v2 E1): the top-level `timing` and `usage`, every seed summary's
    `timing` (results[*].timing), its annotation's lookahead counters (results[*].annotation.direct_reads,
    .keys_mapped) and the caller's seed name with its mismatch flag (results[*].seed.seed_id, .seed_id_mismatch,
    VERIFY-L VL16). Everything else is kept as it is: the strategy (output.timing and annotation.batch_kmers
    included), the capabilities, the release, the rest of the seed summary (validated_seed_id, length_bp, ...)."""
    out = {k: v for k, v in env.items() if k not in ENVELOPE_DROP}
    if isinstance(out.get("results"), list):
        out["results"] = [_summary_view(r) for r in out["results"]]
    return out


def canonical(obj: Any) -> bytes:
    """The suite's canonical JSON (suite_common.canonical_json): sorted keys, no whitespace, UTF-8."""
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")


def envelope_digest(text: Any) -> Optional[str]:
    """sha256 (hex) of the canonical JSON of the J envelope without its per-request values (envelope_view); None for a
    document without J."""
    env = envelope(text)
    if env is None:
        return None
    return hashlib.sha256(canonical(envelope_view(env))).hexdigest()


def digests(text: Any) -> dict:
    """FORMAT-v8 §3's `doc` block of a document: {"body_digest", "envelope_digest", "seed_index", "bytes"}."""
    text = _text(text)
    f = frame(text)
    return {"body_digest": body_digest(text), "envelope_digest": envelope_digest(text),
            "seed_index": f["seed_index"], "bytes": len(text.encode("utf-8"))}
