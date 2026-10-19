#!/usr/bin/env python3
"""The pre-checks of an area G scenario (manifest v8; DESIGN-v8 §2-§4, FORMAT-v8 §5): before the five ways' answers
are compared, whether the ways processed one document, and whether the scenario still is what the manifest froze.

    graphlet_precheck.py RUN_DIR [--out FILE]     (default RUN_DIR/arms/precheck/graphlet_prechecks.jsonl)

One line per scenario with a `graphlet` block, each check an outcome of the suite's vocabulary and what it saw:
  body            the body digests of W1, W2 and W4 (doc.body_digest of each way's records: the sha256 of the body
                  with H's seed_index masked, the J line and the Z record left out) are one; otherwise fail, naming
                  the pair and the first differing record of the documents the ways kept (body_records: the same
                  reading, line by line);
  envelope        the envelope digests (doc.envelope_digest: J without timing and usage) of W1, W2 and W4 are one;
                  otherwise fail, naming the first differing key of the J lines without those values
                  (graphlet_doc.envelope_view; capabilities.feature_level and strategy among the keys);
  twin            checks.graphlet_arms of W1's document (the suite's own MGT reader) against the arms of the
                  `detail: full` twin: the arms, each segment's parents, from_bp, length_bp and sequence (where
                  both state one), each path's segments, length_bp and end_reasons; otherwise fail, naming the first
                  difference;
  strategy        the MCP's echo (traverse_status.overrides: the strategy the service sent) equals
                  request.strategy, and the twin's strategy is request.strategy with output.detail "full"; otherwise
                  fail, naming the drift (the path and both values);
  preconditions   each fact of graphlet.preconditions on W1's document (PRECONDITION_FACTS); otherwise
                  precondition_failed, naming the fact and the value seen;
  library_commit  every W3 call's library.commit (as the service reported it) is the deployed one (what the MCP
                  runner read per invocation, mcp_runs.jsonl / run.json `library_commit`) and the manifest's pin;
                  otherwise fail; not_exposed when the op plan sends no call to W3.
A failed body, envelope or twin pre-check fails the scenario as a determinism finding (DESIGN-v8 §2, §3: with the
rehearsal's headroom a body mismatch is a finding, and a twin that differs leaves W5's answers unattributable); a
failed strategy or library_commit check fails it under its own name. A check that cannot be made (a way recorded
nothing, a document or the twin not kept, no echo recorded) is inconclusive: like every inconclusive it blocks the
gate (the owner, 2026-10-14). The scenario's outcome: fail > precondition_failed > inconclusive > pass.

Where the inputs live (FORMAT-v8 §3 names the records; the rest is lane C's reading, listed as FORMAT questions):
  records         arms/<way dir>/graphlet_ops.jsonl (WAY_DIRS), a later line of one (scenario, op) winning;
  documents       arms/<way dir>/docs/<scenario>.mgt.gz (or .mgt): the text the way processed (W1/W2 the standalone
                  text of the engine's answer); W4 without a copy of its own reads the MCP's own exported file,
                  arms/mcp/docs/<scenario>.mgt.gz (FORMAT v2 C1);
  twin            arms/direct-local/docs/<scenario>.full.json.gz: the engine's `detail: full` answer (the body, or
                  the normalised {http, body} answer);
  echo            the scenario's MCP job record, arms/mcp/graphlet_jobs.jsonl {scenario, job_id, handle, overrides,
                  service_commit, library_commit} (FORMAT v2 C3), else the MCP execution line of the scenario (its raw
                  record's traverse_status call, payload.overrides, as v7's runner records it);
  pin             graphlet.library_commit of the scenario (FORMAT v2 D1, covered by the manifest's content_sha256);
                  the job's library_commit is held to it too.
A failed seed (FORMAT v2 D3: preconditions {"no_graphlet": true}) passes its body, envelope and twin checks when no
way has a document and the twin walked no arm either.
The suite imports nothing of metagraph.traverse here: the digests a record states are graphlet_doc's (lane L's MGT
document reader), recomputed from a kept document to hold the record to it; the record-level reading that names a
difference (body_records, envelope_of, first_difference) is this module's own.
"""

import argparse
import copy
import gzip
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import suite_common as sc  # noqa: E402

try:   # lane L's reader of the MGT documents (body_digest, envelope_digest); the checks work without it
    import graphlet_doc  # noqa: E402
except ImportError:   # pragma: no cover - before lane L's commit lands
    graphlet_doc = None

WAYS = ("W1", "W2", "W3", "W4", "W5")
LIBRARY_WAYS = ("W1", "W2", "W4")
# FORMAT-v8 §3: the directory of each way under arms/ (drift: the report-only lane of the branch head's library)
WAY_DIRS = {"W1": "dl-lib", "W2": "dr-lib", "W3": "mcp", "W4": "mcp-lib", "W5": "oracle", "drift": "drift-lib"}
WAY_NAMES = {"W1": "dl+lib", "W2": "dr+lib", "W3": "mcp", "W4": "mcp+lib", "W5": "oracle", "drift": "drift"}
OPS_FILE = "graphlet_ops.jsonl"
PRECHECK_DIR = "precheck"
PRECHECK_FILE = "graphlet_prechecks.jsonl"
JOBS_FILE = "graphlet_jobs.jsonl"
CHECKS = ("body", "envelope", "twin", "strategy", "preconditions", "library_commit")
# the checks whose failure is a determinism finding (DESIGN-v8 §2, §3)
DETERMINISM_CHECKS = ("body", "envelope", "twin")
# a pre-check's outcomes (the suite's vocabulary; not_exposed: library_commit without a W3 call)
OUTCOMES = ("pass", "fail", "inconclusive", "precondition_failed", "not_exposed")
SCENARIO_ORDER = ("fail", "precondition_failed", "inconclusive", "pass")
# the per-request values of a J line, left out where its keys are named (DESIGN-v8 §2)
PER_REQUEST_KEYS = ("timing", "usage")
# H's field holding the seed's position in the request (SPEC-core §7.5; masked: a batch position, not content)
H_SEED_INDEX_FIELD = 11
TEXT_CUT = 160


# ---------------------------------------------------------------------------------------------------------------
# reading a run's area G records

def graphlet_scenarios(manifest):
    """The scenarios with a `graphlet` block (area G), in manifest order."""
    return [s for s in (manifest or {}).get("scenarios", []) or []
            if isinstance(s, dict) and isinstance(s.get("graphlet"), dict)]


def way_records(run_dir, way):
    """{(scenario, op): record} of one way's graphlet_ops.jsonl (a later line of one key wins, as merge_run keeps a
    resumed execution's last line), and the warnings: unparsable lines, a line of another way."""
    path = os.path.join(run_dir, "arms", WAY_DIRS[way], OPS_FILE)
    rows, bad = sc.read_jsonl(path)
    out, warnings = {}, []
    if bad:
        warnings.append("%s: %d unparsable line(s)" % (os.path.relpath(path, run_dir), bad))
    other = 0
    for r in rows:
        if not isinstance(r, dict):
            continue
        if r.get("way") != way:
            other += 1
            continue
        out[(r.get("scenario"), r.get("op"))] = r
    if other:
        warnings.append("%s: %d line(s) of another way left out" % (os.path.relpath(path, run_dir), other))
    return out, warnings


def all_records(run_dir, ways=WAYS + ("drift",)):
    """({way: {(scenario, op): record}}, warnings, {way: sha256 of its file or None})."""
    recs, warnings, shas = {}, [], {}
    for way in ways:
        recs[way], w = way_records(run_dir, way)
        warnings.extend(w)
        shas[way] = sc.file_sha256(os.path.join(run_dir, "arms", WAY_DIRS[way], OPS_FILE))
    return recs, warnings, shas


def _read_text(path):
    try:
        if path.endswith(".gz"):
            with gzip.open(path, "rb") as fh:
                return fh.read().decode("utf-8")
        with open(path, "rb") as fh:
            return fh.read().decode("utf-8")
    except (OSError, EOFError, UnicodeDecodeError):
        return None


def doc_path(run_dir, way, sid):
    """The document way |way| kept for scenario |sid| (FORMAT v2 C1: arms/<way dir>/docs/<sid>.mgt.gz, or .mgt), or
    None; W4 (and W3) without a copy of its own: the MCP's own exported file, arms/mcp/docs/<sid>.mgt.gz, which W4
    reads."""
    dirs = [WAY_DIRS[way]] + (["mcp"] if way in ("W4", "W3") else [])
    for d in dirs:
        for ext in (".mgt.gz", ".mgt"):
            p = os.path.join(run_dir, "arms", d, "docs", sid + ext)
            if os.path.exists(p):
                return p
    return None


def read_doc(run_dir, way, sid):
    p = doc_path(run_dir, way, sid)
    return _read_text(p) if p else None


def read_twin(run_dir, sid):
    """The `detail: full` twin's answer body (arms/direct-local/docs/<sid>.full.json[.gz]; a normalised answer
    {http, body} gives its body), or (None, why)."""
    for ext in (".full.json.gz", ".full.json"):
        p = os.path.join(run_dir, "arms", "direct-local", "docs", sid + ext)
        if not os.path.exists(p):
            continue
        text = _read_text(p)
        try:
            doc = json.loads(text) if text is not None else None
        except ValueError:
            doc = None
        if not isinstance(doc, dict):
            return None, "the twin %s does not read as a JSON object" % os.path.relpath(p, run_dir)
        if "http" in doc and "body" in doc:
            if doc.get("http") != 200:
                return None, "the twin was answered http %s" % doc.get("http")
            doc = doc.get("body")
        if not isinstance(doc, dict) or not isinstance(doc.get("results"), list) or not doc["results"]:
            return None, "the twin's answer states no results"
        return doc, None
    return None, "no twin kept (arms/direct-local/docs/%s.full.json.gz)" % sid


# ---------------------------------------------------------------------------------------------------------------
# JSON pointers and differences (RFC 6901)

def pointer_escape(key):
    return str(key).replace("~", "~0").replace("/", "~1")


def pointer_get(obj, pointer):
    """(found, value) of an RFC 6901 pointer into |obj| ("" the whole value)."""
    if pointer == "":
        return True, obj
    if not isinstance(pointer, str) or not pointer.startswith("/"):
        return False, None
    cur = obj
    for raw in pointer[1:].split("/"):
        key = raw.replace("~1", "/").replace("~0", "~")
        if isinstance(cur, dict):
            if key not in cur:
                return False, None
            cur = cur[key]
        elif isinstance(cur, list):
            if not key.isdigit() or (len(key) > 1 and key[0] == "0") or int(key) >= len(cur):
                return False, None
            cur = cur[int(key)]
        else:
            return False, None
    return True, cur


def leaf_equal(a, b):
    """Two JSON leaves equal as JSON values (compare.json_leaf_equal: a boolean only the same boolean, numbers by
    value)."""
    if isinstance(a, bool) or isinstance(b, bool):
        return isinstance(a, bool) and isinstance(b, bool) and a == b
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        return a == b
    return type(a) is type(b) and a == b


MISSING = object()


def first_difference(a, b, pointer=""):
    """(pointer, value in a, value in b) of the first difference of two JSON values in sorted-key order, or None;
    a key on one side only gives MISSING on the other."""
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b), key=str):
            p = pointer + "/" + pointer_escape(k)
            if k not in a or k not in b:
                return p, a.get(k, MISSING), b.get(k, MISSING)
            d = first_difference(a[k], b[k], p)
            if d:
                return d
        return None
    if isinstance(a, list) and isinstance(b, list):
        for i, (x, y) in enumerate(zip(a, b)):
            d = first_difference(x, y, "%s/%d" % (pointer, i))
            if d:
                return d
        if len(a) != len(b):
            i = min(len(a), len(b))
            return "%s/%d" % (pointer, i), a[i] if i < len(a) else MISSING, b[i] if i < len(b) else MISSING
        return None
    if isinstance(a, (dict, list)) or isinstance(b, (dict, list)):
        return (pointer, a, b) if not (type(a) is type(b) and a == b) else None
    return None if leaf_equal(a, b) else (pointer, a, b)


def compact(value, limit=TEXT_CUT):
    """A JSON value as one short text ("(absent)" for MISSING)."""
    if value is MISSING:
        return "(absent)"
    text = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return text if len(text) <= limit else text[:limit - 3] + "..."


def difference_dict(diff, a_name="a", b_name="b"):
    """A first_difference as a record's JSON: {pointer, <a_name>, <b_name>} with compact values."""
    if not diff:
        return None
    return {"pointer": diff[0] or "/", a_name: compact(diff[1]), b_name: compact(diff[2])}


def without_keys(obj, keys):
    """|obj| with every key of |keys| left out, at every depth."""
    if isinstance(obj, dict):
        return dict((k, without_keys(v, keys)) for k, v in obj.items() if k not in keys)
    if isinstance(obj, list):
        return [without_keys(v, keys) for v in obj]
    return obj


# ---------------------------------------------------------------------------------------------------------------
# the documents, read in the suite's words

def body_records(text):
    """The records of an MGT document as the body identity reads them: every line but the J line (the second line,
    when it starts with "J ") and the closing Z record, H's seed_index field masked."""
    lines = text.split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    out = []
    for i, line in enumerate(lines):
        if i == 1 and line.startswith("J "):
            continue
        if i == len(lines) - 1 and (line == "Z" or line.startswith("Z ")):
            continue
        if i == 0 and line.startswith("H "):
            f = line.split(" ")
            if len(f) > H_SEED_INDEX_FIELD:
                f[H_SEED_INDEX_FIELD] = "<seed_index>"
            line = " ".join(f)
        out.append(line)
    return out


def first_body_difference(text_a, text_b):
    """{record (1-based, among body_records), type, a, b} of the first record the two documents differ in, or
    None."""
    ra, rb = body_records(text_a), body_records(text_b)
    for i in range(max(len(ra), len(rb))):
        x = ra[i] if i < len(ra) else None
        y = rb[i] if i < len(rb) else None
        if x != y:
            kind = (x or y or "").split(" ", 1)[0]
            cut = lambda s: "(absent: %d records)" % (len(ra) if s is x else len(rb)) if s is None else \
                (s if len(s) <= TEXT_CUT else s[:TEXT_CUT - 3] + "...")  # noqa: E731
            return {"record": i + 1, "type": kind, "a": cut(x), "b": cut(y)}
    return None


def envelope_of(text):
    """The J line of an MGT document as a JSON object, or None (no J line, or not an object)."""
    first = text.find("\n")
    if first < 0 or not text.startswith("J ", first + 1):
        return None
    end = text.find("\n", first + 1)
    try:
        env = json.loads(text[first + 3:end if end >= 0 else len(text)])
    except ValueError:
        return None
    return env if isinstance(env, dict) else None


def lane_l_digests(text):
    """(body_digest, envelope_digest) of a document by graphlet_doc (lane L's reader), or None without it."""
    if graphlet_doc is None or text is None:
        return None
    try:
        return graphlet_doc.body_digest(text), graphlet_doc.envelope_digest(text)
    except Exception:   # a document lane L's reader cannot read: no recomputed digests
        return None


# ---------------------------------------------------------------------------------------------------------------
# the pre-checks

def _check(outcome, detail, **extra):
    return dict({"outcome": outcome, "detail": detail}, **extra)


def scenario_ops(scenario):
    g = scenario.get("graphlet") or {}
    return [op for op in g.get("ops") or [] if isinstance(op, dict) and op.get("id")]


def op_ways(scenario, op):
    """The ways an op runs in: its own `ways`, else the scenario's ways.default, else all five."""
    g = scenario.get("graphlet") or {}
    ways = op.get("ways") if isinstance(op.get("ways"), list) else ((g.get("ways") or {}).get("default")
                                                                    if isinstance(g.get("ways"), dict) else None)
    return [w for w in WAYS if w in (ways if isinstance(ways, list) else WAYS)]


def _way_digests(recs, sid, way, key):
    """The distinct doc.<key> values way |way|'s records of |sid| state (None left out), in order."""
    out = []
    for (s, _), r in sorted(recs.get(way, {}).items(), key=lambda kv: str(kv[0][1])):
        if s != sid:
            continue
        d = (r.get("doc") or {}).get(key) if isinstance(r.get("doc"), dict) else None
        if d and d not in out:
            out.append(d)
    return out


def _has_records(recs, sid, way):
    return any(s == sid for (s, _) in recs.get(way, {}))


def digest_check(run_dir, recs, sid, key, docs, differ):
    """The body or envelope identity across W1, W2 and W4 (|key| body_digest / envelope_digest; |differ|(text_a,
    text_b) names the first difference of two kept documents)."""
    digests, missing, several = {}, [], []
    for way in LIBRARY_WAYS:
        got = _way_digests(recs, sid, way, key)
        if not got:
            missing.append(way)
        elif len(got) > 1:
            several.append("%s states %d different %ss" % (way, len(got), key))
        digests[way] = got[0] if got else None
    # the documents kept, against the digests their records state (graphlet_doc where it loads)
    recomputed = {}
    for way in LIBRARY_WAYS:
        d = lane_l_digests(docs.get(way))
        if d is not None:
            recomputed[way] = d[0] if key == "body_digest" else d[1]
    stale = ["%s's records state %s %s, its kept document's is %s" % (w, key, str(digests[w])[:12], str(v)[:12])
             for w, v in sorted(recomputed.items()) if digests.get(w) and v != digests[w]]
    if several or stale:
        return _check("fail", "; ".join(several + stale), digests=digests)
    present = [w for w in LIBRARY_WAYS if digests.get(w)]
    if len(set(digests[w] for w in present)) > 1:
        ref = present[0]
        other = next(w for w in present if digests[w] != digests[ref])
        pair = "%s-%s" % (ref, other)
        diff = None
        if docs.get(ref) is not None and docs.get(other) is not None:
            diff = differ(docs[ref], docs[other])
        detail = "%s: %s %s vs %s" % (pair, key, digests[ref][:12], digests[other][:12])
        if diff:
            detail += "; first difference %s" % (
                "record %d (%s): %s | %s" % (diff["record"], diff["type"], diff["a"], diff["b"]) if "record" in diff
                else "%s: %s vs %s" % (diff["pointer"], diff[ref], diff[other]) if ref in diff
                else json.dumps(diff, sort_keys=True))
        else:
            detail += "; the documents are not both kept (arms/<way>/docs/%s.mgt.gz): the record is not named" % sid
        return _check("fail", detail, digests=digests, pair=pair, difference=diff)
    if missing:
        return _check("inconclusive", "no %s recorded by %s" % (key, ", ".join(missing)), digests=digests)
    return _check("pass", "%s %s on W1, W2, W4" % (key, digests["W1"][:12]), digests=digests)


def body_check(run_dir, recs, sid, docs):
    return digest_check(run_dir, recs, sid, "body_digest", docs, first_body_difference)


def envelope_view(env):
    """A J envelope without its per-request values: graphlet_doc.envelope_view (the definition envelope_digest
    hashes: the top-level timing and usage, every results[i].timing) where lane L's reader loads, else every `timing`
    and `usage` key at any depth."""
    if graphlet_doc is not None and hasattr(graphlet_doc, "envelope_view"):
        return graphlet_doc.envelope_view(env)
    return without_keys(env, PER_REQUEST_KEYS)


def _envelope_differ(ref, other):
    def differ(text_a, text_b):
        ea, eb = envelope_of(text_a), envelope_of(text_b)
        if ea is None or eb is None:
            return {"pointer": "/", ref: "(no J line)" if ea is None else "(J)",
                    other: "(no J line)" if eb is None else "(J)"}
        return difference_dict(first_difference(envelope_view(ea), envelope_view(eb)), ref, other) or \
            {"pointer": "/", ref: "(equal without the per-request values)", other: "(digests differ)"}
    return differ


def envelope_check(run_dir, recs, sid, docs):
    bare = [w for w in LIBRARY_WAYS if docs.get(w) is not None and envelope_of(docs[w]) is None]
    if bare:
        # DESIGN-v8 §2 (D-R1): the envelope travels with the body; without it summary and compare differ by construction
        return _check("fail", "%s's document carries no J line (the engine's envelope)" % ", ".join(bare))
    digests = dict((w, (_way_digests(recs, sid, w, "envelope_digest") or [None])[0]) for w in LIBRARY_WAYS)
    present = [w for w in LIBRARY_WAYS if digests.get(w)]
    ref = present[0] if present else "W1"
    other = next((w for w in present if digests[w] != digests[ref]), "W4")
    return digest_check(run_dir, recs, sid, "envelope_digest", docs, _envelope_differ(ref, other))


def _arm_shape(arm):
    segs = [dict((k, s.get(k)) for k in ("parents", "from_bp", "length_bp")) for s in arm.get("segments") or []]
    seqs = [s.get("sequence") for s in arm.get("segments") or []]
    paths = [{"segments": p.get("segments"), "length_bp": p.get("length_bp"),
              "end_reasons": p.get("end_reasons") or {}} for p in arm.get("paths") or []]
    return segs, seqs, paths


def twin_difference(graphlet_text, full_body, checks_mod):
    """The first difference between the arms of an MGT document (checks.graphlet_arms) and the `detail: full`
    twin's (results[0].arms): {what, arm, ..., graphlet, full}, or None."""
    g_arms = checks_mod.graphlet_arms(graphlet_text)
    f_arms = ((full_body.get("results") or [{}])[0] or {}).get("arms") or {}
    if sorted(g_arms) != sorted(f_arms):
        return {"what": "arms", "graphlet": sorted(g_arms), "full": sorted(f_arms)}
    for side in sorted(f_arms):
        gs, gq, gp = _arm_shape(g_arms[side])
        fs, fq, fp = _arm_shape(f_arms[side] or {})
        for what, a, b in (("segments", gs, fs), ("paths", gp, fp)):
            d = first_difference(a, b)
            if d:
                return {"what": what, "arm": side, "pointer": d[0] or "/", "graphlet": compact(d[1]),
                        "full": compact(d[2])}
        for i, (x, y) in enumerate(zip(gq, fq)):
            if x is not None and y is not None and x != y:
                return {"what": "sequence", "arm": side, "pointer": "/segments/%d/sequence" % i,
                        "graphlet": compact(x), "full": compact(y)}
    return None


def twin_check(run_dir, sid, docs, checks_mod):
    text = docs.get("W1")
    if text is None:
        return _check("inconclusive", "W1 kept no document (arms/%s/docs/%s.mgt.gz)" % (WAY_DIRS["W1"], sid))
    full, why = read_twin(run_dir, sid)
    if full is None:
        return _check("inconclusive", why)
    if checks_mod is None or not hasattr(checks_mod, "graphlet_arms"):
        return _check("inconclusive", "the suite's MGT reader (checks.graphlet_arms) did not load")
    try:
        diff = twin_difference(text, full, checks_mod)
    except Exception as e:   # a body the reader rejects: the twin cannot be compared with it
        return _check("fail", "W1's document does not read (%s: %s)" % (type(e).__name__, str(e)[:200]))
    if diff:
        where = "arms" if diff["what"] == "arms" else "%s arm, %s %s" % (diff["arm"], diff["what"], diff["pointer"])
        return _check("fail", "the graphlet differs from its full twin: %s: graphlet %s, full %s" % (
            where, compact(diff["graphlet"]) if diff["what"] == "arms" else diff["graphlet"],
            compact(diff["full"]) if diff["what"] == "arms" else diff["full"]), difference=diff)
    return _check("pass", "W1's graphlet has the full twin's arms, segments, paths, lengths and end reasons")


def _payload_value(v):
    """A raw payload value as the runner keeps it: a JSON object, or its JSON text."""
    if isinstance(v, str):
        try:
            return json.loads(v)
        except ValueError:
            return v
    return v


# FORMAT v2 C3: the MCP job record, one line per scenario in arms/mcp/graphlet_jobs.jsonl: job id, handle (the
# record's UUID), the traverse_status.overrides echo, service_commit, library_commit
ECHO_KEYS = ("overrides", "overrides_echo", "echo")


def mcp_job(run_dir, sid):
    """The scenario's MCP job record (the last line naming it), or {}."""
    rows, _ = sc.read_jsonl(os.path.join(run_dir, "arms", "mcp", JOBS_FILE))
    out = {}
    for r in rows:
        if isinstance(r, dict) and r.get("scenario") == sid:
            out = r
    return out


def mcp_echo(run_dir, sid):
    """(the strategy the service sent for scenario |sid| (traverse_status.overrides), where it was read) or (None,
    why): arms/mcp/graphlet_jobs.jsonl first, else the scenario's first MCP execution's raw record."""
    echo = None
    job = mcp_job(run_dir, sid)
    for key in ECHO_KEYS:
        if job.get(key) is not None:
            echo = _payload_value(job[key])
            break
    if isinstance(echo, dict):
        return echo, "arms/mcp/%s" % JOBS_FILE
    arm_dir = os.path.join(run_dir, "arms", "mcp")
    lines, _ = sc.read_jsonl(os.path.join(arm_dir, "executions.jsonl"))
    mine = sorted((r for r in lines if isinstance(r, dict) and r.get("scenario") == sid),
                  key=lambda r: r.get("repeat") if isinstance(r.get("repeat"), int) else 99)
    for line in mine:
        doc = sc.read_raw_record(line, arm_dir=arm_dir)
        for call in reversed((doc or {}).get("calls") or []):
            if isinstance(call, dict) and call.get("tool") == "traverse_status":
                ov = _payload_value((call.get("payload") or {}).get("overrides"))
                if isinstance(ov, dict):
                    return ov, "the raw record of %s mcp r%s" % (sid, line.get("repeat"))
    return None, "no traverse_status echo recorded for %s (arms/mcp/%s, or the MCP execution's raw calls)" % (
        sid, JOBS_FILE)


def _dotted(pointer):
    return ".".join(p.replace("~1", "/").replace("~0", "~") for p in pointer.split("/")[1:]) or "<root>"


def strategy_check(run_dir, scenario):
    sid = scenario.get("id")
    req = scenario.get("request") if isinstance(scenario.get("request"), dict) else {}
    strategy = req.get("strategy")
    if not isinstance(strategy, dict):
        return _check("inconclusive", "the scenario's request states no strategy")
    twin = ((scenario.get("graphlet") or {}).get("twin") or {}).get("strategy")
    drifts = []
    if isinstance(twin, dict):
        want = copy.deepcopy(strategy)
        want.setdefault("output", {})["detail"] = "full"
        d = first_difference(twin, want)
        if d:
            drifts.append({"what": "twin", "path": _dotted(d[0]), "twin": compact(d[1]), "request": compact(d[2])})
    echo, where = mcp_echo(run_dir, sid)
    if echo is None:
        if drifts:
            return _check("fail", "the twin's strategy drifts from request.strategy at %s (twin %s, request %s); %s" % (
                drifts[0]["path"], drifts[0]["twin"], drifts[0]["request"], where), drifts=drifts)
        return _check("inconclusive", where)
    d = first_difference(echo, strategy)
    if d:
        drifts.insert(0, {"what": "echo", "path": _dotted(d[0]), "echo": compact(d[1]), "request": compact(d[2])})
    if drifts:
        text = "; ".join("strategy drift at %s: %s %s, request %s" % (
            x["path"], "the MCP echo" if x["what"] == "echo" else "the twin",
            x.get("echo", x.get("twin")), x["request"]) for x in drifts)
        return _check("fail", text, drifts=drifts, echo_source=where)
    return _check("pass", "the MCP echo (%s) is request.strategy%s" % (
        where, ", the twin's strategy is it with output.detail full" if isinstance(twin, dict) else ""),
        echo_source=where)


# graphlet.preconditions: fact -> (what it reads, how it is judged); every other key is a manifest problem
PRECONDITION_FACTS = {
    "arms": "the arms the document carries, exactly these",
    "walks_min": "{arm: n}: at least n walks (leaves) on the arm",
    "segments_min": "n (over the arms) or {arm: n}: at least n segments",
    "merges_min": "n (over the arms) or {arm: n}: at least n segments with more than one parent",
    "labels_present": "[name]: each name in the label dictionary (L records)",
    "walk_bp_min": "{arm: bp}: the arm's longest walk at least bp long",
    "bases": "true / false: the segments carry their bases (a sequences: false retrieval does not)",
    "no_graphlet": "true: a failed seed, W1 has no graphlet (FORMAT v2 D3)",
}


def _per_arm_or_total(want, values, name, out):
    if isinstance(want, dict):
        for arm, n in sorted(want.items()):
            seen = values.get(arm, 0)
            out.append({"fact": "%s.%s" % (name, arm), "want": n, "seen": seen,
                        "outcome": "pass" if isinstance(n, int) and seen >= n else "fail"})
    else:
        seen = sum(values.values())
        out.append({"fact": name, "want": want, "seen": seen,
                    "outcome": "pass" if isinstance(want, int) and seen >= want else "fail"})


def precondition_facts(text, pre, checks_mod):
    """[{fact, want, seen, outcome pass | fail | inconclusive}] of graphlet.preconditions on an MGT document."""
    out = []
    arms = checks_mod.graphlet_arms(text)
    for key, want in sorted(pre.items()):
        if key == "arms":
            seen = sorted(arms)
            out.append({"fact": key, "want": want, "seen": seen,
                        "outcome": "pass" if isinstance(want, list) and sorted(want) == seen else "fail"})
        elif key == "walks_min":
            _per_arm_or_total(want, dict((a, len(arms[a]["paths"])) for a in arms), key, out)
        elif key == "segments_min":
            _per_arm_or_total(want, dict((a, len(arms[a]["segments"])) for a in arms), key, out)
        elif key == "merges_min":
            _per_arm_or_total(want, dict((a, sum(1 for s in arms[a]["segments"] if len(s["parents"]) > 1))
                                         for a in arms), key, out)
        elif key == "labels_present":
            names = set(x.get("name") for x in checks_mod.graphlet_label_dict(text))
            for name in want if isinstance(want, list) else [want]:
                out.append({"fact": "labels_present", "want": name, "seen": name in names,
                            "outcome": "pass" if name in names else "fail"})
        elif key == "walk_bp_min":
            for arm, bp in sorted((want or {}).items()) if isinstance(want, dict) else []:
                seen = max([p["length_bp"] for p in (arms.get(arm) or {}).get("paths") or []] or [0])
                out.append({"fact": "walk_bp_min.%s" % arm, "want": bp, "seen": seen,
                            "outcome": "pass" if isinstance(bp, int) and seen >= bp else "fail"})
        elif key == "bases":
            segs = [s for a in arms.values() for s in a["segments"]]
            seen = bool(segs) and all("sequence" in s for s in segs)
            out.append({"fact": key, "want": want, "seen": seen, "outcome": "pass" if seen is want else "fail"})
        else:
            out.append({"fact": key, "want": want, "seen": None, "outcome": "inconclusive",
                        "detail": "not a precondition this suite reads (%s)" % ", ".join(sorted(PRECONDITION_FACTS))})
    return out


def no_graphlet(scenario):
    """Whether the scenario is a failed seed (FORMAT v2 D3: graphlet.preconditions {"no_graphlet": true})."""
    pre = (scenario.get("graphlet") or {}).get("preconditions")
    return isinstance(pre, dict) and pre.get("no_graphlet") is True


def preconditions_check(scenario, docs, checks_mod):
    pre = (scenario.get("graphlet") or {}).get("preconditions")
    if not pre:
        return _check("pass", "the scenario states no preconditions", facts=[])
    if not isinstance(pre, dict):
        return _check("inconclusive", "graphlet.preconditions is not an object")
    if no_graphlet(scenario):
        # a failed seed: the fact is that W1 has no graphlet (no document; its calls answer no_graphlet)
        fact = {"fact": "no_graphlet", "want": True, "seen": docs.get("W1") is None,
                "outcome": "pass" if docs.get("W1") is None else "fail"}
        if fact["outcome"] == "fail":
            return _check("precondition_failed", "no_graphlet: want true, seen a document of W1 (the seed has a "
                                                 "graphlet now)", facts=[fact])
        return _check("pass", "the seed has no graphlet on W1, as the scenario states", facts=[fact])
    text = docs.get("W1")
    if text is None:
        return _check("inconclusive", "W1 kept no document to check the preconditions on")
    if checks_mod is None:
        return _check("inconclusive", "the suite's MGT reader (checks.graphlet_arms) did not load")
    try:
        facts = precondition_facts(text, pre, checks_mod)
    except Exception as e:
        return _check("inconclusive", "W1's document does not read (%s: %s)" % (type(e).__name__, str(e)[:200]))
    failed = [f for f in facts if f["outcome"] == "fail"]
    unread = [f for f in facts if f["outcome"] == "inconclusive"]
    if failed:
        return _check("precondition_failed", "; ".join("%s: want %s, seen %s" % (
            f["fact"], compact(f["want"]), compact(f["seen"])) for f in failed), facts=facts)
    if unread:
        return _check("inconclusive", "; ".join("%s: %s" % (f["fact"], f["detail"]) for f in unread), facts=facts)
    return _check("pass", "%d fact(s) hold on W1's document" % len(facts), facts=facts)


def library_pin(scenario, manifest=None):
    """(the manifest's library_commit pin, where it is stated) for a scenario, or (None, None): FORMAT v2 D1,
    graphlet.library_commit per scenario, inside content_sha256 (a top-level block would be outside it)."""
    pin = (scenario.get("graphlet") or {}).get("library_commit")
    if isinstance(pin, str) and pin.strip():
        return pin.strip(), "graphlet.library_commit"
    return None, None


def commit_matches(recorded, expected):
    """A recorded commit is the expected one: equal ignoring case, or one an abbreviation (a prefix of at least 7 hex
    digits) of the other (gate.commit_matches)."""
    a, b = str(recorded or "").strip().lower(), str(expected or "").strip().lower()
    if len(a) < 7 or len(b) < 7:
        return False
    return a.startswith(b) or b.startswith(a)


def deployed_library_commits(run_dir):
    """{commit: invocations} the MCP runner recorded as the deployed library (library_commit per invocation in
    arms/mcp/mcp_runs.jsonl, else run.json invocations[], else run.json), and the number of invocations without one."""
    lines, _ = sc.read_jsonl(os.path.join(run_dir, "arms", "mcp", "mcp_runs.jsonl"))
    records = [r for r in lines if isinstance(r, dict)]
    if not records:
        run = sc.read_json(os.path.join(run_dir, "arms", "mcp", "run.json")) or {}
        records = [r for r in run.get("invocations") or [] if isinstance(r, dict)] or ([run] if run else [])
    commits, without = {}, 0
    for r in records:
        c = r.get("library_commit")
        if isinstance(c, str) and c.strip():
            commits[c.strip()] = commits.get(c.strip(), 0) + 1
        else:
            without += 1
    return commits, without


def library_commit_check(run_dir, recs, scenario, manifest, deployed=None):
    sid = scenario.get("id")
    w3_ops = [op["id"] for op in scenario_ops(scenario) if "W3" in op_ways(scenario, op)]
    if not w3_ops:
        return _check("not_exposed", "the op plan sends no call to W3")
    pin, where = library_pin(scenario, manifest)
    commits, without = deployed if deployed is not None else deployed_library_commits(run_dir)
    seen, missing, bad = {}, [], []
    for op in w3_ops:
        r = recs.get("W3", {}).get((sid, op))
        if r is None:
            continue
        c = (r.get("library") or {}).get("commit") if isinstance(r.get("library"), dict) else None
        if not c:
            missing.append(op)
            continue
        seen[c] = seen.get(c, 0) + 1
        if pin and not commit_matches(c, pin):
            bad.append("%s: %s" % (op, c[:12]))
        elif any(not commit_matches(c, d) for d in commits):
            bad.append("%s: %s, deployed %s" % (op, c[:12], ", ".join(d[:12] for d in sorted(commits))))
    extra = {"pin": pin, "pin_source": where, "deployed": commits, "calls": seen}
    if bad:
        return _check("fail", "W3 call(s) report another library_commit than the %s: %s" % (
            "manifest's pin %s" % pin[:12] if pin else "deployed one", "; ".join(bad[:10])), **extra)
    if pin and any(not commit_matches(d, pin) for d in commits):
        return _check("fail", "the deployed library_commit %s is not the manifest's pin %s" % (
            ", ".join(d[:12] for d in sorted(commits)), pin[:12]), **extra)
    job_commit = mcp_job(run_dir, sid).get("library_commit")
    extra["job"] = job_commit
    if pin and job_commit and not commit_matches(job_commit, pin):
        return _check("fail", "the scenario's MCP job reports the library_commit %s, not the manifest's pin %s" % (
            str(job_commit)[:12], pin[:12]), **extra)
    if missing:
        return _check("inconclusive", "%d W3 call(s) report no library_commit (%s)" % (
            len(missing), ", ".join(missing[:10])), **extra)
    if not seen:
        return _check("inconclusive", "no W3 record of the scenario's %d W3 op(s)" % len(w3_ops), **extra)
    if not pin:
        return _check("inconclusive", "the manifest pins no library_commit (graphlet.library_commit)", **extra)
    if not commits:
        return _check("inconclusive", "the MCP runner recorded no deployed library_commit (mcp_runs.jsonl)", **extra)
    return _check("pass", "every W3 call (%d) reports the pinned and deployed library_commit %s" % (
        sum(seen.values()), pin[:12]), **extra)


def scenario_outcome(checks):
    """(outcome, finding): fail > precondition_failed > inconclusive > pass over the checks; finding "determinism"
    when a body, envelope or twin check failed, "strategy drift" / "library_commit" for those."""
    outs = dict((k, (checks.get(k) or {}).get("outcome")) for k in CHECKS)
    if "fail" in outs.values():
        if any(outs[k] == "fail" for k in DETERMINISM_CHECKS):
            finding = "determinism"
        elif outs["strategy"] == "fail":
            finding = "strategy drift"
        else:
            finding = "library_commit"
        return "fail", finding
    for o in ("precondition_failed", "inconclusive"):
        if o in outs.values():
            return o, None
    return "pass", None


def failed_seed_checks(run_dir, recs, sid, docs):
    """The body, envelope and twin checks of a failed seed (FORMAT v2 D3): no way has a document or states a digest,
    and the twin found no graphlet either (its result states no arms)."""
    kept = [w for w in LIBRARY_WAYS if docs.get(w) is not None or _way_digests(recs, sid, w, "body_digest")]
    ident = _check("fail", "a failed seed, yet %s has a document" % ", ".join(kept)) if kept else \
        _check("pass", "a failed seed: no way has a document")
    full, why = read_twin(run_dir, sid)
    if full is None:
        twin = _check("pass", "a failed seed: the twin has no graphlet either (%s)" % why) \
            if why and "states no results" in why else _check("inconclusive", why)
    elif any((r or {}).get("arms") for r in full.get("results") or []):
        twin = _check("fail", "a failed seed, yet the full twin walked arms")
    else:
        twin = _check("pass", "a failed seed: the twin's result has no arms either")
    return ident, dict(ident), twin


def precheck_scenario(run_dir, scenario, recs, manifest, checks_mod, deployed=None):
    """The pre-check record of one area G scenario (FORMAT-v8 §5)."""
    sid = scenario.get("id")
    docs = dict((w, read_doc(run_dir, w, sid)) for w in LIBRARY_WAYS)
    if no_graphlet(scenario) and docs.get("W1") is None:
        body, envelope, twin = failed_seed_checks(run_dir, recs, sid, docs)
    else:
        body, envelope = body_check(run_dir, recs, sid, docs), envelope_check(run_dir, recs, sid, docs)
        twin = twin_check(run_dir, sid, docs, checks_mod)
    checks = {
        "body": body,
        "envelope": envelope,
        "twin": twin,
        "strategy": strategy_check(run_dir, scenario),
        "preconditions": preconditions_check(scenario, docs, checks_mod),
        "library_commit": library_commit_check(run_dir, recs, scenario, manifest, deployed),
    }
    outcome, finding = scenario_outcome(checks)
    rec = {"scenario": sid, "outcome": outcome, "finding": finding,
           "documents": dict((w, os.path.relpath(doc_path(run_dir, w, sid), run_dir) if doc_path(run_dir, w, sid)
                              else None) for w in LIBRARY_WAYS)}
    rec.update(checks)
    return rec


def load_checks(run_dir=None):
    """The suite's checks module (its MGT reader), as the comparison loads it (compare.judging_checks: the run's
    snapshot where it has one), else the one next to this script; None when none loads."""
    try:
        import compare as cmp
        mod = cmp.load_checks_module(cmp.judging_checks(run_dir)["dir"]) if run_dir else None
        if mod is not None and hasattr(mod, "graphlet_arms"):
            return mod
    except Exception:
        pass
    try:
        import checks as mod
        return mod
    except Exception:
        return None


def prechecks(run_dir, manifest=None, recs=None, checks_mod=None):
    """[pre-check record] of every area G scenario of the run (manifest: RUN_DIR/manifest.json by default)."""
    manifest = manifest if manifest is not None else sc.read_json(os.path.join(run_dir, "manifest.json"))
    if recs is None:
        recs, _, _ = all_records(run_dir)
    checks_mod = checks_mod if checks_mod is not None else load_checks(run_dir)
    deployed = deployed_library_commits(run_dir)
    return [precheck_scenario(run_dir, s, recs, manifest, checks_mod, deployed)
            for s in graphlet_scenarios(manifest)]


def precheck_path(run_dir):
    return os.path.join(run_dir, "arms", PRECHECK_DIR, PRECHECK_FILE)


def write_prechecks(run_dir, records, path=None):
    path = path or precheck_path(run_dir)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        for r in records:
            fh.write(json.dumps(r, sort_keys=True, separators=(",", ":"), ensure_ascii=False) + "\n")
    os.replace(tmp, path)
    return path


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("run_dir")
    ap.add_argument("--out", default=None, help="default RUN_DIR/arms/precheck/graphlet_prechecks.jsonl")
    a = ap.parse_args(argv)
    manifest = sc.read_json(os.path.join(a.run_dir, "manifest.json"))
    if not isinstance(manifest, dict):
        sys.stderr.write("graphlet_precheck.py: %s has no manifest.json\n" % a.run_dir)
        return 2
    records = prechecks(a.run_dir, manifest)
    path = write_prechecks(a.run_dir, records, a.out)
    counts = {}
    for r in records:
        counts[r["outcome"]] = counts.get(r["outcome"], 0) + 1
    sys.stdout.write("graphlet pre-checks: %d scenario(s) %s; wrote %s\n" % (
        len(records), ", ".join("%s %d" % kv for kv in sorted(counts.items())) or "none", path))
    return 0 if all(r["outcome"] == "pass" for r in records) else 1


if __name__ == "__main__":
    sys.exit(main())
