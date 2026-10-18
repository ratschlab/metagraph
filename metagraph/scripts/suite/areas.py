"""The scenario catalogue of DESIGN §3, one generator per area (R, T, P, L, X, M, Q, C, D, N; S is the async
session's). Each scenario names its ground truth, the exact host body, the MCP mapping (manifest_lib), its
expectations (check names of checks.py with their arguments), a timing class and tags.

Expectations are what the CURRENT deployed code answers per the SPECs and the C++ source (texts quoted from
metagraph/src and the stored fixtures); a scenario that depends on a host capability is built from the build's
capabilities snapshot (manifest_lib.HOST_MISMATCH_MAX: the protein substitution cap), one whose expectation the smoke
run has not confirmed yet and that rests on reading the code is tagged code-derived.
"""

from __future__ import annotations

import copy
import json
import os
from typing import Any, Dict, List, Optional, Tuple

import seqtools as st
from manifest_lib import (ABSENT_LABEL, HOST_INDEX, HOST_MISMATCH_MAX, SERVER_MAX_SEED_LABELS, allowed_statuses,
                          SERVICE_DISCOVER_MAX_LABELS, Ctx, DB, K,
                          PATTERN_MAX_CONTEXTS, RESOLVE_MAX_QUERY_BP, acgt_flank, trie_heads, trie_index,
                          PATTERN_TIME_BUDGET_MS, RESOLVE_TIME_BUDGET_MS, WALK_PRESET, pattern_request, relation,
                          _pattern_kinds, resolve_request, scenario, source_contexts, traverse_request,
                          walk_expectations)
from truth import TRUTH_DIR

# A search the host may not finish within its time budget on a loaded machine (the README's machine-dependent) asks
# for a 200 too: checks.deadline_allowed allows exactly the SPEC's 503 `deadline` in its place (SPEC-pattern-search
# §6 / §7.6, SPEC-labeled-traversal-core §4.5 with a request budget; A-02 of the review of lane S5 and A2-02 of the
# review of lane S5b: a status of [200, 503] passed any 503, a proxy page or `loading` too; no builder states it).
OK = {"check": "status", "args": {"http": 200}}
INV = {"check": "invariants", "args": {}}
CC = {"check": "counts_consistent", "args": {}}

# genomes with few relatives in RefSeq 97 (their windows carry few columns: labels "all" and discovery stay small)
LOW_CARRIER = ["buc", "mge", "mpn", "ser", "mja", "sso", "tth", "pfa", "cal", "hvo"]
SMALL_WHOLE = ["lam", "t4", "hsv1", "hiv1", "hsam", "drem", "dmem", "scem"]


def ex(check: str, **args) -> dict:
    return {"check": check, "args": args}


def empirical(t: dict, claim: str) -> dict:
    """The truth |t| with the part of the scenario's claim that the records cannot certify (a graph-wide or a
    column-wide absence, which other records of the index decide), named with the index it was observed on: a failure
    on another index asks to re-validate the expectation before it reads as a regression. The scenario carries the
    tag index-empirical."""
    return dict(t, index_empirical={"index_fp": HOST_INDEX.get("index_fp"), "claim": claim})


LABEL_DICT_NOT_VIA_MCP = "the MCP traverse answers carry no label dictionary"


def mcp_not_implemented_check(sc: dict, check: str, reason: str) -> None:
    """State in the manifest that the MCP arm cannot answer one property check of |sc| yet
    (mcp.not_implemented_checks [{check, reason}]): the scenario's other checks stay tested via MCP, this one stays
    required on the direct arms (gate_policy derives the MCP exemption, the report lists it as not implemented yet via
    MCP, never as untested)."""
    if "tool" in sc["mcp"]:
        sc["mcp"].setdefault("not_implemented_checks", []).append({"check": check, "reason": reason})


def mcp_by_design_check(sc: dict, check: str, case: str) -> None:
    """State in the manifest that the service answers one property check of |sc| otherwise by design, the table's
    |case| (mcp.by_design_checks [{check, case}]): required on the direct arms, service-side by design on MCP."""
    sc["mcp"].setdefault("by_design_checks", []).append({"check": check, "case": case})


def err400(text: Optional[str] = None, code: Optional[str] = None, prefix: Optional[str] = None,
           contains: Optional[List[str]] = None, http: int = 400) -> dict:
    a: Dict = {"http": http}
    if code is not None:
        a["code"] = code
    if text is not None:
        a["text"] = text
    if prefix is not None:
        a["text_prefix"] = prefix
    if contains:
        a["text_contains"] = contains
    return {"check": "status", "args": a}


def _ids(prefix: str):
    n = [0]

    def nxt() -> str:
        n[0] += 1
        return "%s-%04d" % (prefix, n[0])
    return nxt


def _record_runs(rec, q: str) -> List[List[int]]:
    """The runs ([start, end) of k-mer indices) of the query k-mers |rec|'s fetched sequence holds on its forward strand:
    on a whole record, exactly what its header label supports on this basic graph."""
    out: List[List[int]] = []
    for i in range(len(q) - K + 1):
        if q[i:i + K] in rec.seq:
            if out and out[-1][1] == i:
                out[-1][1] = i + 1
            else:
                out.append([i, i + 1])
    return out


def _distant(ctx: Ctx, rec, i: int):
    others = [r for r in ctx.positives if r.kingdom != rec.kingdom]
    return others[i % len(others)]


def _columns_sample() -> List[str]:
    with open(os.path.join(TRUTH_DIR, "index", "columns_first1200.txt")) as f:
        return [l.strip() for l in f if l.strip() and not l.startswith("#")]


# The warm-ups' ids, fixed per target (a warm-up added later takes the next free number; no id moves)
WARM_IDS = {"T-repeat-0001": "R-warm-0001", "T-repeat-0002": "R-warm-0002", "T-repeat-0003": "R-warm-0003",
            "T-repeat-0004": "R-warm-0004", "T-labels-0003": "R-warm-0005", "N-hard-0022": "R-warm-0006",
            "T-walk-0009": "R-warm-0007"}
# A walk whose warm-up is split in two, each resolve with its own 60,000 ms budget: the seed with its left extent
# (WARM_IDS), then the seed's last k-1 bases with its right extent (this id). T-repeat-0002: its seed rows alone took
# 17.6 s cold in v3 and more than 22 s in v4, one right-arm row 11.7 s, so the whole 295 bp could exceed one budget
# under v3's load (VERIFY-V5A's re-verification)
SPLIT_WARM_IDS = {"T-repeat-0002": "R-warm-0008"}
# what a walk's warm-up reads beyond the walk's exact minimum on each arm (the rows of the steps that decide it)
WARM_MARGIN_BP = 45


def warm_up(ctx: Ctx, sid: str, query: str, label: str, truth: dict, target: str, why: str,
            tags: Tuple[str, ...] = ()) -> dict:
    """A warm-up of scenario |target|: |query| resolved against |label| (the record's own label: every k-mer
    supported, a positive check of its own under the default 60,000 ms budget), which reads the annotation rows the
    target reads (a header label the coordinate tuples, a column its rows), placed right before the target
    (build_manifest.place_warm_ups), so that the target's first execution, the MCP arm's cold one, meets them warm:
    v4's cold first calls under the Logan load timed out where the warm repeats took a second (triage-v4 G24)."""
    assert len(query) <= RESOLVE_MAX_QUERY_BP, "a warm-up the MCP can send (at most %d bp)" % RESOLVE_MAX_QUERY_BP
    t = dict(truth, warm_up=why)
    return scenario(sid, "R", "warm-up of %s: %s bp resolved against %s, every k-mer supported" % (
                        target, format(len(query), ","), label), "POST /resolve",
                    resolve_request(query, labels=[label]), t,
                    [OK, INV, ex("resolve_runs", label=label, full=True)], tags=("warm-up", "positive") + tuple(tags),
                    warms=target)


def region_warm_up(ctx: Ctx, rec, t: dict, left: int, right: int, label: str, target: str, why: str,
                   sid: Optional[str] = None) -> dict:
    """warm_up (id |sid|, else WARM_IDS[target]) of the record region a walk from the seed |t| must read: the seed,
    |left| bases before it and |right| after it (cut at a base outside ACGT), against the walk's label."""
    s0 = t["start"] - rec.offset
    a, _, _ = acgt_flank(rec, s0, t["length"], left)
    _, b, _ = acgt_flank(rec, s0, t["length"], right)
    region = ctx.truth_at(rec, a, b - a, flank=0)
    region["flank"] = ""
    return warm_up(ctx, sid or WARM_IDS[target], region["window"], label, region, target, why,
                   tags=("header-label",) if label == rec.accession else ("column-label",))


def minimum_warm_up(ctx: Ctx, sc: dict, label: str, why: str) -> List[dict]:
    """region_warm_up of what a header walk must read to decide walk_in_record: the seed and, on each requested arm,
    its exact minimum (truth.walk[side].min_bp, _annotate_walks) plus WARM_MARGIN_BP, within the radius; past the
    minimum a time stop no longer hides the property. A cold resolve of a multi-copy region costs seconds per row
    (VERIFY-V5A: T-repeat-0002's 71 seed rows took 17.6 s cold, one row 11.7 s), so the warm-up reads no more. A
    target of SPLIT_WARM_IDS gets two: the seed with its left extent, then the seed's last k-1 bases (the right arm's
    first k-mers start there) with its right extent."""
    t = sc["truth"]
    rec = ctx.rec(t["genome"])
    radius = sc["request"]["strategy"]["bounds"]["max_extension_bp"]
    walk = t.get("walk") or {}
    ext = dict((side, min(radius, (walk.get(side) or {}).get("min_bp", radius) + WARM_MARGIN_BP)) if side in walk else
               (side, 0) for side in ("left", "right"))
    if sc["id"] not in SPLIT_WARM_IDS:
        return [region_warm_up(ctx, rec, t, ext["left"], ext["right"], label, sc["id"], why)]
    tail = dict(t, start=t["start"] + t["length"] - (K - 1), length=K - 1)
    return [region_warm_up(ctx, rec, t, ext["left"], 0, label, sc["id"], why + " (the seed and its left arm)"),
            region_warm_up(ctx, rec, tail, 0, ext["right"], label, sc["id"],
                           why + " (the seed's last %d bases and its right arm)" % (K - 1), SPLIT_WARM_IDS[sc["id"]])]


# ================================================================================================= R: /resolve
def gen_R(ctx: Ctx) -> List[dict]:
    out: List[dict] = []
    rng = ctx.rng("R")
    recs = ctx.positives
    nid = _ids("R-pos")
    for li, L in enumerate([31, 32, 50, 100, 200, 1000, 2000]):
        for j in range(2):
            rec = recs[(li * 2 + j * 7) % len(recs)]
            t = ctx.window(rng, rec, L)
            for kind in ("column", "header"):
                lab = rec.column if kind == "column" else rec.accession
                req = resolve_request(t["window"], labels=[lab])
                exp = [OK, INV, ex("resolve_runs", label=lab, full=True),
                       ex("resolve_graph_runs", runs=[[0, L - K + 1]])]
                out.append(scenario(nid(), "R", "%d bp of %s against its %s label %s" % (L, rec.key, kind, lab),
                                    "POST /resolve", req, t, exp, tags=["positive", kind + "-label"]
                                    + (["smoke"] if (L == 100 and j == 0 and kind == "column") else [])))
    nid = _ids("R-multi")
    for i in range(6):
        rec = recs[(i * 5 + 3) % len(recs)]
        other = _distant(ctx, rec, i)
        t = ctx.window(rng, rec, 100)
        labels = [rec.column, other.column, rec.accession]
        req = resolve_request(t["window"], labels=labels)
        exp = [OK, INV, ex("resolve_runs", label=rec.column, full=True),
               ex("resolve_runs", label=rec.accession, full=True), ex("label_not_full", label=other.column)]
        t["note"] = ("label_not_full assumes the unrelated taxon (another kingdom) lacks one of this random, "
                     "DUST-filtered 100 bp window's k-mers; a window of a gene conserved across kingdoms would break it")
        out.append(scenario(nid(), "R", "%s window against its taxid, its header and the unrelated %s" % (rec.key, other.key),
                            "POST /resolve", req, t, exp, tags=["positive", "multi-label"]))
    nid = _ids("R-junction")
    for i in range(4):
        ra = recs[(i * 3 + 1) % len(recs)]
        rb = _distant(ctx, ra, i + 2)
        ta = ctx.window(rng, ra, 100)
        tb = ctx.window(rng, rb, 100)
        q = ta["window"] + tb["window"]
        req = resolve_request(q, labels=[ra.column, rb.column])
        # the exact split: A's k-mers run past 70 while B's first bases repeat A's own continuation, B's start before
        # 100 while A's last bases repeat what precedes B in its record (computed from the records)
        a_next = ra.seq[ta["start"] - ra.offset + 100:ta["start"] - ra.offset + 100 + K]
        b_prev = rb.seq[max(0, tb["start"] - rb.offset - K):tb["start"] - rb.offset]
        cp = next((x for x in range(min(len(a_next), 100)) if tb["window"][x] != a_next[x]), min(len(a_next), 100))
        cs = next((x for x in range(min(len(b_prev), 100)) if ta["window"][-1 - x] != b_prev[-1 - x]),
                  min(len(b_prev), 100))
        a_end, b_start = 70 + min(cp, K - 1), 100 - min(cs, K - 1)
        n = len(q) - K + 1
        # The runs are exact for the records' headers (the k-mers each whole record holds, which is all a header
        # supports on this basic graph); a taxid column holds the column's other records too, which may carry the
        # junction's bases (a variant of the column, P5-v4-content-03: pae's column holds 1,762 records), so a
        # column's one run from the query's end reaches into the junction within it, never past it: the k-mers of
        # the other genome are in no record of the column (another kingdom's; empirical, the index's answers).
        req = resolve_request(q, labels=[ra.column, rb.column, ra.accession, rb.accession])
        exp = [OK, INV, ex("resolve_runs", label=ra.accession, runs=_record_runs(ra, q)),
               ex("resolve_runs", label=rb.accession, runs=_record_runs(rb, q)),
               ex("resolve_runs", label=ra.column, run_end_in=[a_end, 100], max_kmers=100),
               ex("resolve_runs", label=rb.column, run_start_in=[70, b_start], max_kmers=n - 70)]
        truth = empirical({"a": ta, "b": tb, "junction_bp": 100, "a_run_end": a_end, "b_run_start": b_start},
                          "the taxid columns hold no k-mer of the other genome (their other records, another "
                          "kingdom's); the headers' runs are computed from the whole records")
        out.append(scenario(nid(), "R", "junction of %s and %s: A's run ends at %d, B's starts at %d (the headers' "
                            "exactly, from the records; the columns' within the junction)" % (
                                ra.key, rb.key, a_end, b_start),
                            "POST /resolve", req, truth, exp, tags=["positive", "junction", "index-empirical"]))
    # SPEC-labeled-traversal-core §4.2 ranks a discovery by supported k-mers (ties by column id, then seq_id) and
    # truncates only above max_labels, so the source, which supports every in-graph k-mer, is left out of a truncated
    # list only in a tie at full support (resolve_discover_contains reads that from labels_truncated, lane S7c). At
    # limit 1 the v4 windows of buc (column) and ser (header) were left out in such ties (8 and 2 labels at full
    # support), which tested the tie and never the listing: their limit is the v4 tie plus one, so the source must be
    # listed. Every discovery names max_labels and kind, at most 100: the MCP sends its limit and kind as they are and
    # refuses above 100, the host's own default is 1,000 (the async session's constraint for v5).
    full_support_tie_v4 = {"buc": {"total": 9, "full_support": 8}, "ser": {"total": 7, "full_support": 2}}
    nid = _ids("R-discover")
    for i, key in enumerate(["buc", "mge", "mpn", "ser", "mja", "sso"]):
        rec = ctx.rec(key)
        t = ctx.window(rng, rec, 100)
        kind = "column" if i % 2 == 0 else "header"
        limit = [1, 20, 100][i % 3]
        tie = full_support_tie_v4.get(key)
        if tie:
            limit = tie["full_support"] + 1
            t = empirical(dict(t, discover_v4=dict(tie, max_labels=1)),
                          "the labels carrying the window at full support number %d (index 907ec986's answers in v3 "
                          "and v4): the limit lies above them" % tie["full_support"])
        req = resolve_request(t["window"], discover={"max_labels": limit, "kind": kind})
        lab = rec.column if kind == "column" else rec.accession
        exp = [OK, INV, ex("resolve_discover_contains", label=lab)]
        out.append(scenario(nid(), "R", "discover %s labels (limit %d%s) on a %s window: the source listed (ranked by "
                            "supported k-mers; a truncated list leaves it out only in a tie at full support)" % (
                                kind, limit, ", above v4's full-support tie of %d" % tie["full_support"] if tie else "",
                                key),
                            "POST /resolve", req, t, exp, timing_class="slow",
                            tags=["discover"] + (["smoke"] if i == 0 else []) + (["index-empirical"] if tie else [])))
    nid = _ids("R-trace")
    for i, key in enumerate(["eco", "sau", "lam"]):
        rec = ctx.rec(key)
        t = ctx.window(rng, rec, 200)
        req = resolve_request(t["window"], labels=[rec.accession], support="trace")
        exp = [OK, INV, ex("resolve_runs", label=rec.accession, full=True)]
        out.append(scenario(nid(), "R", "trace support of %s's header on its own 200 bp: one contiguous run" % key,
                            "POST /resolve", req, t, exp, tags=["trace", "header-label"]))
    # A stated stop (SPEC-labeled-traversal-core §4.5) is a 200 only when the row batch in flight at the work deadline
    # (budget - 250 ms, the fixed reserve) ends, and the answer is built, before the budget; past it the answer is the
    # 503. Cold or loaded rows make every batch long: no answer of v3 or v4 stated a stop (triage-v4 V05). So each
    # budget scenario follows its complete twin in manifest order (R-twin-0002 / -0004 under the 60,000 ms budget,
    # its warm-up) and its budget puts the work deadline inside the warm work, after the first batches (64 rows, then
    # doubling up to 4,096): pae's 10 kb took 343 ms warm (v4 R-twin-0004 direct-remote), so 420 ms (the deadline at
    # 170 ms falls in the 4,096-row batch, which ends at about 276 ms at that rate); eco's 10 kb never completed on
    # this host (v3 and v4: 21-45 s), its second half took 1.37 s warm (R-twin-0003), so the query is that half under
    # 480 ms (the owner, 2026-10-15; VERIFY-V5A's re-verification): the work deadline at 230 ms falls in a batch that
    # ends before the budget for warm totals of about 0.28 to 2.48 s (the batch model's window, stop_window_warm_ms),
    # where 700 ms gave a stop only for 1.13 to 1.75 s, provided the unpolled k-mer lookup ends within the 230 ms (pae's
    # whole 10 kb, lookup included, took 343 ms warm). pae's window is about 0.21 to 21.8 s. Whether a stop happens
    # still depends on the batch times, i.e. on the load and on what is warm; a complete answer reads inconclusive (no
    # stop, nothing tested). The scenarios stay in the deploy gate, required on direct-local (the owner, 2026-10-15),
    # one repeat each: a repeat decides on its own, so a second would double what must stop (the builder keeps them off
    # the paired plan's warm-repeat index, build_manifest.keep_off_warm_repeat).
    nid = _ids("R-budget")
    for i, (key, budget, half, warm_ms, window) in enumerate((("eco", 480, True, 1370, [284, 2480]),
                                                              ("pae", 420, False, 343, [209, 21800]))):
        rec = ctx.rec(key)
        t = ctx.window(rng, rec, 10000, unique=False, max_dust=5.0)
        if half:
            t = ctx.truth_at(rec, t["start"] - rec.offset + 5000, 5000)
        L = t["length"]
        req = resolve_request(t["window"], labels=[rec.column], time_budget_ms=budget)
        t2 = dict(t, flank="", stop_design={
            "rule": "the work deadline (budget - 250 ms) falls inside the warm work, in a row batch that ends before "
                    "the budget at the warm rate measured",
            "warm_work_ms": warm_ms, "work_deadline_ms": budget - 250,
            "stop_window_warm_ms": window,
            "window_assumes": "the k-mer lookup, which the deadline does not poll, ends within the work deadline; "
                              "rows in batches of 64 doubling up to 4,096 at a uniform warm rate"})
        exp = [OK, INV, ex("resolve_runs", label=rec.column, full=True, prefix_if_stopped=True)]
        out.append(scenario(nid(), "R", "%s bp of %s under a %d ms budget after its complete twin: a stated stop answers "
                            "the prefix, the record's own label supporting every k-mer of it" % (format(L, ","), key, budget),
                            "POST /resolve", req, t2, exp, timing_class="slow", tags=["budget", "machine-dependent"]))
    # R-edge
    nid = _ids("R-edge")
    eco = ctx.rec("eco")
    t = ctx.window(rng, eco, 31)
    out.append(scenario(nid(), "R", "exactly k bases: one k-mer, supported", "POST /resolve",
                        resolve_request(t["window"], labels=[eco.column]), t,
                        [OK, INV, ex("resolve_runs", label=eco.column, full=True), ex("resolve_graph_runs", runs=[[0, 1]])],
                        tags=["edge"]))
    t = ctx.window(rng, eco, 30, unique=False)
    out.append(scenario(nid(), "R", "k-1 bases: refused", "POST /resolve",
                        resolve_request(t["window"], labels=[eco.column]), t,
                        [err400("Query shorter than k")], tags=["edge", "refusal", "smoke"]))
    t = ctx.window(rng, eco, 10001, unique=False, max_dust=5.0)
    t["flank"] = ""
    out.append(scenario(nid(), "R", "10,001 bp: above the server's max_query_bp (10,000 on staging)", "POST /resolve",
                        resolve_request(t["window"], labels=[eco.column]), t,
                        [err400("request.sequence: longer than the server limit of 10000 bp")], tags=["edge", "refusal"]))
    t = ctx.window(rng, eco, 200)
    out.append(scenario(nid(), "R", "200 bp above the request's bounds.max_query_bp 100", "POST /resolve",
                        resolve_request(t["window"], labels=[eco.column], max_query_bp=0,
                                        bounds={"max_query_bp": 100, "time_budget_ms": RESOLVE_TIME_BUDGET_MS}), t,
                        [err400("request.sequence: longer than bounds.max_query_bp")], tags=["edge", "refusal"]))
    for code, label in (("N", "N"), ("R", "IUPAC R")):
        t = ctx.window(rng, eco, 100)
        q = t["window"][:50] + code + t["window"][51:]
        out.append(scenario(nid(), "R", "an %s at position 50: the k-mers holding it are not in the graph" % label,
                            "POST /resolve", resolve_request(q, labels=[eco.column]), t,
                            [OK, INV, ex("resolve_graph_runs", runs=[[0, 20], [51, 70]]),
                             ex("resolve_runs", label=eco.column, kmers=39)],
                            tags=["edge", "code-derived"]))
    t = ctx.window(rng, eco, 100)
    out.append(scenario(nid(), "R", "lower case: resolved as upper case", "POST /resolve",
                        resolve_request(t["window"].lower(), labels=[eco.column]), t,
                        [OK, INV, ex("resolve_runs", label=eco.column, full=True)], tags=["edge", "code-derived"]))
    t = ctx.window(rng, eco, 100)
    out.append(scenario(nid(), "R", "an unknown label: refused naming it", "POST /resolve",
                        resolve_request(t["window"], labels=[eco.column, ABSENT_LABEL]), t,
                        [err400("Label not found in the annotation: '%s'" % ABSENT_LABEL)], tags=["edge", "refusal"]))
    out.append(scenario(nid(), "R", "a duplicate label: refused (LabelQuery), not merged", "POST /resolve",
                        resolve_request(t["window"], labels=[eco.column, eco.column]), t,
                        [err400("Duplicate label '%s'" % eco.column)], tags=["edge", "refusal"]))
    cols = [c for c in _columns_sample() if c != eco.column][:1000] + [eco.column]
    out.append(scenario(nid(), "R", "1,001 explicit labels: the host has no cap (the service's is 1,000)", "POST /resolve",
                        resolve_request(t["window"], labels=cols), t,
                        [OK, INV, ex("resolve_runs", label=eco.column, full=True)], timing_class="slow", tags=["edge"]))
    out.append(scenario(nid(), "R", "an empty sequence: refused", "POST /resolve",
                        resolve_request("", labels=[eco.column]), None, [err400("request.sequence is required")],
                        tags=["edge", "refusal"]))
    out.append(scenario(nid(), "R", "discover together with labels: refused", "POST /resolve",
                        {"sequence": t["window"], "labels": [eco.column], "discover": {"max_labels": 5}}, t,
                        [err400("request: specify exactly one of 'labels' or 'discover'")], tags=["edge", "refusal"]))
    out.append(scenario(nid(), "R", "neither labels nor discover: refused", "POST /resolve",
                        {"sequence": t["window"], "bounds": {"time_budget_ms": RESOLVE_TIME_BUDGET_MS}}, t,
                        [err400("request: specify exactly one of 'labels' or 'discover'")], tags=["edge", "refusal"]))
    q = st.random_dna(ctx.rng("R", "absent"), 200)
    out.append(scenario(nid(), "R", "a random 200 bp: no k-mer in the graph", "POST /resolve",
                        resolve_request(q, labels=[eco.column]), {"random": True},
                        [OK, INV, ex("resolve_absent")], tags=["negative"]))
    out.append(scenario(nid(), "R", "an unknown field: refused naming it", "POST /resolve",
                        dict(resolve_request(t["window"], labels=[eco.column]), bogus=1), t,
                        [err400("request: unknown field 'bogus'")], tags=["edge", "refusal"]))
    out.append(scenario(nid(), "R", "graph on a single-graph server: refused", "POST /resolve",
                        dict(resolve_request(t["window"], labels=[eco.column]), graph="x"), t,
                        [err400("Bad request: this server hosts a single graph; remove the 'graph' / 'graph_path' field")],
                        tags=["edge", "refusal"]))
    out.append(scenario(nid(), "R", "graphs on a single-graph server: refused", "POST /resolve",
                        dict(resolve_request(t["window"], labels=[eco.column]), graphs=["x"]), t,
                        [err400("Bad request: this server hosts a single graph; remove the 'graphs' field")],
                        tags=["edge", "refusal"]))
    out.append(scenario(nid(), "R", "a body that is not an object: refused", "POST /resolve", [t["window"]], None,
                        [err400("request: expected an object")], tags=["edge", "refusal", "protocol"]))
    # an object cut before its closing brace: no JSON reader accepts it (/resolve's reader takes a trailing comma or a
    # comment, jsoncpp's defaults: G12; only /pattern promises RFC 8259 strictness)
    out.append(scenario(nid(), "R", "a body that is not JSON (an object cut before its end): refused", "POST /resolve",
                        None, None, [err400(prefix="Bad json received: ")], tags=["edge", "refusal", "protocol"],
                        body_raw='{"sequence": "%s", "labels": ["%s"]' % (t["window"], eco.column)))
    out.append(scenario(nid(), "R", "time_budget_ms not above the finalisation reserve: refused", "POST /resolve",
                        resolve_request(t["window"], labels=[eco.column], time_budget_ms=250), t,
                        [err400(prefix="request.bounds.time_budget_ms: expected more than the finalisation reserve of 250")],
                        tags=["edge", "refusal"]))
    out.append(scenario(nid(), "R", "time_budget_ms above the server cap: lowered and stated", "POST /resolve",
                        resolve_request(t["window"], labels=[eco.column], time_budget_ms=70000), t,
                        [OK, INV, ex("clamped", field="bounds.time_budget_ms", requested=70000, effective=60000),
                         ex("resolve_runs", label=eco.column, full=True)], tags=["edge", "clamp"]))
    out.append(scenario(nid(), "R", "discover.max_labels 0: refused naming the field", "POST /resolve",
                        resolve_request(t["window"][:60], discover={"max_labels": 0, "kind": "column"}), t,
                        [err400(prefix="request.discover.max_labels: out of range [1, ")], tags=["edge", "refusal"]))
    # the service's discovery cap is 1,000 since the v5 start (it was 100): 150 is answered on every arm; the source,
    # eco's column, is listed or ties the kept labels at full support (v4: 1,443 of 2,753 columns at 30 of 30 k-mers)
    out.append(scenario(nid(), "R", "discover.max_labels 150: answered on every arm, the source listed or in a tie at "
                        "full support", "POST /resolve",
                        resolve_request(t["window"][:60], discover={"max_labels": 150, "kind": "column"}),
                        t, [OK, INV, ex("resolve_discover_contains", label=eco.column)],
                        tags=["edge", "discover", "machine-dependent"]))
    out.append(scenario(nid(), "R", "time_budget_ms 400 (above the reserve, below the service's 500): the host answers",
                        "POST /resolve", resolve_request(t["window"], labels=[eco.column], time_budget_ms=400), t,
                        [OK, INV, ex("budget_respected")], tags=["edge", "split-case", "budget", "machine-dependent"]))
    # above the service's cap: the MCP refuses it (invalid_arguments), the host, which has no cap, answers (the table's
    # resolve.discover_limit_above_cap_host_answers)
    out.append(scenario(nid(), "R", "discover.max_labels 1,001 (above the service's cap of 1,000): the host answers, "
                        "the MCP refuses", "POST /resolve",
                        resolve_request(t["window"][:60], discover={"max_labels": SERVICE_DISCOVER_MAX_LABELS + 1,
                                                                     "kind": "column"}),
                        t, [OK, INV], tags=["edge", "split-case", "discover", "machine-dependent"]))
    # the reverse complement of a window none of whose reverse k-mers its record holds (Ctx.rc_negative_window: the
    # whole record, k-mer by k-mer): the record's header supports none of the query's k-mers (on this basic graph a
    # k-mer and its reverse complement are different nodes). The check discriminates only where the reverse k-mers are
    # in the graph through other records (an assembly in the other orientation): genomes with many RefSeq records.
    # R-rc-0001's S. aureus window was in the graph whole in v3 and v4; v4's R-rc-0002 (an M. jannaschii window) was
    # absent from the graph, which made it a plain absence: E. coli takes its place
    nid = _ids("R-rc")
    for key, evidence in (("sau", "this window's reverse complement: graph_runs [[0, 70]] on every arm of v3 and v4"),
                          ("eco", "E. coli's many RefSeq records (6,181 headers carried T-labels-0004's window in v4); "
                                  "the reverse complement of an eco 25-mer met 66 forward contexts (v4 N-usecase-0006)")):
        rec = ctx.rec(key)
        t = ctx.rc_negative_window(ctx.rng("R", "rc", key), rec, 100)
        out.append(scenario(nid(), "R", "the reverse complement of a window of %s none of whose reverse k-mers its "
                            "record holds: its header supports none of it" % key,
                            "POST /resolve", resolve_request(st.revcomp(t["window"]), labels=[rec.accession]),
                            dict(t, strand="-", reverse_kmers_in_record=0, reverse_in_graph_evidence=evidence),
                            [OK, INV, ex("resolve_runs", label=rec.accession, kmers=0)],
                            tags=["negative", "reverse-complement", "header-label"]))
    # the label order: the same labels in another order answer each label alike (resolve_runs reads them by name)
    nid = _ids("R-order")
    rec = ctx.rec("hvo")
    other = _distant(ctx, rec, 1)
    t = ctx.window(ctx.rng("R", "order"), rec, 100)
    out.append(scenario(nid(), "R", "labels in reverse order (header, unrelated, taxid): each answered as alone",
                        "POST /resolve", resolve_request(t["window"], labels=[rec.accession, other.column, rec.column]), t,
                        [OK, INV, ex("resolve_runs", label=rec.column, full=True),
                         ex("resolve_runs", label=rec.accession, full=True)], tags=["positive", "multi-label"]))
    return out


# ================================================================================================= T: /traverse
def _tseed(seq: str, labels: Optional[List[str]] = None) -> dict:
    s = {"sequence": seq}
    if labels is not None:
        s["labels"] = labels
    return s


def gen_T(ctx: Ctx) -> List[dict]:
    out: List[dict] = []
    rng = ctx.rng("T")
    recs = [r for r in ctx.positives]
    # The radii are the class's (100, 1,000 or 10,000 bp) except where the Logan load showed the walk cannot finish
    # within its budget (triage-v4 V06, G24; item 5 of lane V5, as the owner amended it on 2026-10-15). The windows
    # are drawn with the class's margin, so that every window stays as built.
    # - The three 1,000 bp header walks that a cold first call left inconclusive in v4 (T-walk-0004, -0014, -0022:
    #   the header's coordinates read at 36-58 bp/s, and the MCP's walk preset holds the budget at 15,000 ms) walk
    #   300 bp, at most 9 s with the seed phase; the other nine keep 1,000 bp (the owner accepts that one may go
    #   inconclusive under load).
    # - The four 10,000 bp column walks that never completed in 15 s in v3 or v4 are split. The gate's scenario keeps
    #   the id, asks for the host's cap, 60,000 ms (the MCP refuses above its preset's 15,000 ms, by design: a direct-
    #   arm walk, whose direct-local repeat 0 is then the cold call), and walks what a cold first call reached in about
    #   30 s (VERIFY-V5A: the frontiers of v3's and v4's successive 15 s calls chained, since v4's direct-local calls
    #   ran after the MCP's had warmed the walk's start): 1,000 / 500 / 5,000 / 1,000 bp. Its long twin
    #   (T-longwalk-nnnn, right after it) is the v4 test: 10,000 bp, under the host's cap of 60,000 ms instead of v4's
    #   15,000, which stopped all 32 executions of v3 and v4 by time (the MCP refuses the cap above its preset, by
    #   design: a direct-arm benchmark); whether it completes depends on the load, so it is budget_stress, left out of
    #   the deploy gate (gate_policy's rule) and kept in the full suite. -0024 and -0030 end at stated branches within
    #   their first 4 kb and passed on every arm in 15 s.
    warm: List[dict] = []                          # warm-ups built in place
    warm_specs: List[Tuple[str, str, str, str]] = []   # warm-ups of a walk's minimum, built after _annotate_walks
    SHORT = {6: 1000, 12: 500, 18: 5000, 36: 1000}
    LONG_TWIN = {6: 10000, 12: 10000, 18: 10000, 36: 10000}
    # walks whose seed row alone took most of the budget cold: a warm-up of the seed and its radius right before
    COLD_ROW_WALKS = {"T-walk-0009": "the seed's row in lambda's column, whose cold read took 15.7 s in v3's first "
                                     "MCP call (time-stopped at 1 bp; VERIFY-V5A)"}
    # the frontier v3's and v4's four successive 15 s calls of the 10 kb walk reached (60 s of cold reads, chained)
    COLD_FRONTIER_60S = {6: (2731, 6398), 12: (1431, 1643), 18: (8147, 8748), 36: (2244, 2255)}
    HEADER_300 = (4, 14, 22)
    nid = _ids("T-walk")
    nid_long = _ids("T-longwalk")
    for i in range(36):
        L = [31, 50, 100, 200][i % 4]
        direction = ["right", "left", "both"][i % 3]
        drawn = 10000 if i % 6 == 5 else [100, 1000][i % 2]
        ext = SHORT.get(i + 1, 300 if i + 1 in HEADER_300 else drawn)
        kind = "column" if i % 3 == 2 else "header"
        rec = recs[(i * 7 + 2) % len(recs)]
        margin = drawn + 200
        if len(rec.seq) < 2 * margin + L:
            rec = ctx.rec("eco")
        tw = ctx.window(rng, rec, L, margin=margin)
        t = ctx.flank_for(rec, tw, ext + 50)
        lab = rec.accession if kind == "header" else rec.column
        exp = [OK, INV, ex("traverse_seed_outcome", walks=["complete", "partial"]),
               ex("walk_in_record", label_scope="record" if kind == "header" else "column")]
        sid = nid()
        if i + 1 not in SHORT:
            req = traverse_request([_tseed(t["window"], [lab])],
                                   {"direction": direction, "bounds": {"max_extension_bp": ext}})
            out.append(scenario(sid, "T", "%d bp seed of %s, %s arm(s), %d bp, constrained to its %s label" % (
                                    L, rec.key, direction, ext, kind),
                                "POST /traverse", req, t, exp, timing_class="slow" if drawn >= 10000 else "fast",
                                tags=["walk", kind + "-label"] + (["smoke"] if i == 1 else [])))
            if sid in COLD_ROW_WALKS:
                warm.append(region_warm_up(ctx, rec, t, ext, ext, lab, sid, COLD_ROW_WALKS[sid]))
            continue
        long_id = nid_long()
        req = traverse_request([_tseed(t["window"], [lab])],
                               {"direction": direction, "bounds": {"max_extension_bp": ext, "time_budget_ms": 60000}})
        out.append(scenario(sid, "T", "%d bp seed of %s, %s arm(s), %s bp, constrained to its %s label, under 60,000 "
                            "ms (the gate's walk; %s walks %s bp)" % (
                                L, rec.key, direction, format(ext, ","), kind, long_id, format(LONG_TWIN[i + 1], ",")),
                            "POST /traverse", req, t, exp, timing_class="slow", tags=["walk", kind + "-label"]))
        big = LONG_TWIN[i + 1]
        front = COLD_FRONTIER_60S[i + 1]
        tl = dict(ctx.flank_for(rec, tw, big + 50), budget_stress={
            "reason": "a %s bp column-label walk under the load: 60 s of cold reads (v3's and v4's four successive "
                      "15 s calls, chained) reached %s / %s bp, so whether one 60,000 ms call completes depends on the "
                      "load: no budget makes it decisive; %s, its window walked %s bp, is the gate's" % (
                          format(big, ","), format(front[0], ","), format(front[1], ","), sid, format(ext, ",")),
            "short_twin": sid, "cold_frontier_60s_v3_v4": list(front)})
        req = traverse_request([_tseed(t["window"], [lab])],
                               {"direction": direction, "bounds": {"max_extension_bp": big, "time_budget_ms": 60000}})
        out.append(scenario(long_id, "T", "%d bp seed of %s, %s arm(s), %s bp, constrained to its %s label, under "
                            "60,000 ms (the long twin of %s; budget_stress)" % (
                                L, rec.key, direction, format(big, ","), kind, sid),
                            "POST /traverse", req, tl, exp, timing_class="slow",
                            tags=["walk", kind + "-label", "budget_stress"]))
    nid = _ids("T-derived")
    for i, key in enumerate(LOW_CARRIER[:6]):
        rec = ctx.rec(key)
        t = ctx.flank_for(rec, ctx.window(rng, rec, 100, margin=1300), 1050)
        req = traverse_request([_tseed(t["window"])], {"bounds": {"max_extension_bp": 1000}})
        exp = [OK, INV, ex("traverse_seed_outcome", walks=["complete", "partial"]),
               ex("walk_in_record", label_scope="derived")]
        out.append(scenario(nid(), "T", "%s seed with its derived carrier set (permit per_hit)" % key, "POST /traverse",
                            req, t, exp, tags=["walk", "derived-labels"]))
    nid = _ids("T-end")
    # records without terminal repeats (HSV-1's and HIV-1's ends are repeated inside them); a 31 bp seed (one k-mer,
    # a validation of few rows) 40-100 bp from the end, so that the walk reaches the end unless a stated stop says
    # otherwise, under the walk preset's 15,000 ms: the MCP refuses a budget above its preset row (v4's 60,000 ms made
    # it refuse all six by design, so the record-end stop was never walked via MCP; P4-builder-truth-01), and the
    # direct walks took 12 ms to 8.9 s
    for i, key in enumerate(["lam", "t4", "hsam", "drem", "dmem", "scem"]):
        rec = ctx.rec(key)
        side = "right" if i % 2 == 0 else "left"
        n = len(rec.seq)
        where = (n - 100 - 31, n - 40 - 31) if side == "right" else (40, 100)
        t = None
        for dust in (2.0, 3.0):
            try:
                t = ctx.window(rng, rec, 31, margin=0, where=where, unique=True, max_dust=dust)
                break
            except ValueError:
                continue
        if t is None:
            ctx.skip("T", "T-end %s %s" % (key, side), "no unique 31 bp window 40-100 bp from the %s end" % side)
            continue
        t = ctx.flank_for(rec, t, 300)
        req = traverse_request([_tseed(t["window"], [rec.accession])],
                               {"direction": side, "bounds": {"max_extension_bp": 200}})
        exp = [OK, INV, ex("walk_end", arm=side, reasons_in=["label_lost", "dead_end"]),
               ex("walk_in_record", label_scope="record")]
        out.append(scenario(nid(), "T", "31 bp seed near the %s end of %s: the walk stops at the record's end" % (side, key),
                            "POST /traverse", req, t, exp, tags=["walk", "record-end"] + (["smoke"] if i == 0 else [])))
    # A repeat's seed and walk read the rows of multi-copy k-mers (rRNA operons, IS elements), which cold under the
    # load took up to 22 s for the seed phase alone (T-repeat-0002's MCP r0 met its attempt's bound, T-repeat-0004's
    # walks ran at 20-53 bp/s; triage-v4 G24): each walk follows a warm-up resolve of its record region under its
    # header, and walks 400 bp, past every branch the record states for these seeds (by 198 bp at the latest); the
    # windows are chosen as before (copies diverging within 950 bp)
    T_REPEAT_RADIUS = 400
    nid = _ids("T-repeat")
    for i, (key, t, copies) in enumerate(_repeat_windows(ctx, rng, n=4, ext=1000)):
        rec = ctx.rec(key)
        R = T_REPEAT_RADIUS
        t = ctx.flank_for(rec, t, R + 50)
        t["copies"] = [{"start": c, "flank": ctx.flank_for(rec, ctx.truth_at(rec, c - rec.offset, t["length"]), R + 50)["flank"],
                        "flank_start": ctx.flank_for(rec, ctx.truth_at(rec, c - rec.offset, t["length"]), R + 50)["flank_start"]}
                       for c in copies if c != t["start"]]
        req = traverse_request([_tseed(t["window"], [rec.accession])], {"bounds": {"max_extension_bp": R}})
        exp = [OK, INV, ex("walk_in_record", label_scope="record", any_copy=True, expect_branching=True)]
        sid = nid()
        out.append(scenario(sid, "T", "seed in a %d-copy repeat of %s: a branch stated where the record branches "
                            "(copies diverge at %s bp), a path follows a copy" % (
                                len(copies), key, "/".join(str(t["copy_divergence_bp"][x]) for x in ("left", "right"))),
                            "POST /traverse", req, t, exp, tags=["walk", "repeat"]))
        warm_specs.append(("minimum", sid, rec.accession,
                           "the coordinate tuples of a multi-copy repeat that the seed and the walk's minimum read"))
    nid = _ids("T-labels")
    eco = ctx.rec("eco")
    for i in range(2):
        rec = ctx.rec(LOW_CARRIER[i])
        t = ctx.flank_for(rec, ctx.window(rng, rec, 100, margin=400), 250)
        req = traverse_request([_tseed(t["window"])], {"bounds": {"max_extension_bp": 200}},
                               extra={"labels": {"mode": "annotate"}, "branching": {"max_label_branches": "unlimited"}})
        out.append(scenario(nid(), "T", "annotate mode on %s: every structural successor, a path follows the record" % rec.key,
                            "POST /traverse", req, t,
                            [OK, INV, ex("walk_in_record", label_scope="annotate")], tags=["labels", "annotate"]))
    # seed_label_kind header vs column on one seed of a multi-record column (the verification's A3: on a
    # single-record column both kinds carry the same k-mers): the derived set holds the record's own label of the
    # kind asked, and only labels of that kind
    # A derived set keeps the first labels.max_seed_labels carriers in (column, seq_id) order (SPEC §6.1). The offline
    # truth cannot count a window's header carriers (it knows the genomes' own records, not which RefSeq records hold
    # a window), so the counts come from the v4 answers (triage-v4 V02): the eco window's header set held 6,181
    # headers (every E. coli assembly; NC_000913.3 falls past the default 1,000, G14), the cal window's 1. Only eco's
    # header kind asks for more labels (the server's maximum, 10,000, the service's too), and with its 6,181 live
    # labels the modelled memory of detail full outgrows the walk preset's 512 MB (a stop at 58 bp in v4), so it asks
    # 4,096 MB (the host states no memory cap: max_memory_mb 0). Both knobs are on the MCP's advanced path or above
    # its preset row: that MCP arm stays not implemented (traverse.advanced_path_disabled); cal's pair walks via MCP.
    # eco's derived sets are large (868 columns, 6,181 headers): its column-kind walk ran at 6-12 bp/s cold under the
    # load (T-labels-0003, v4), so eco's pair walks 50 bp and the column kind follows a warm-up resolve of its region
    # under eco's column (the rows its walk reads); cal's pair walks 200 bp as before
    seed_label_counts = {("eco", "header"): 6181, ("cal", "header"): 1}
    for key in ("eco", "cal"):
        rec = ctx.rec(key)
        t = ctx.flank_for(rec, ctx.window(rng, rec, 100, margin=400), 250)
        radius = 50 if key == "eco" else 200
        for kind in ("column", "header"):
            labels = {"seed_label_kind": kind}
            bounds: Dict[str, Any] = {"max_extension_bp": radius}
            many = seed_label_counts.get((key, kind), 0) > 1000
            tt = t
            if (key, kind) in seed_label_counts:
                tt = dict(t, derived_labels_v4=seed_label_counts[(key, kind)])
            if many:
                labels["max_seed_labels"] = SERVER_MAX_SEED_LABELS
                bounds["max_memory_mb"] = 4096
            req = traverse_request([_tseed(t["window"])], {"bounds": bounds, "labels": labels})
            name = rec.column if kind == "column" else rec.accession
            sid = nid()
            out.append(scenario(sid, "T", "derived %s labels on %s (%d records in its column%s): %s among them" % (
                                    kind, rec.key, rec.records_in_column,
                                    "; max_seed_labels %d, 4,096 MB" % SERVER_MAX_SEED_LABELS if many else "", name),
                                "POST /traverse", req, tt,
                                [OK, INV, ex("walk_in_record", label_scope="derived"),
                                 ex("label_dict_contains", name=name, kind=kind, only_kind=True)],
                                tags=["labels", "derived-labels", "label-kind-pair"]))
            # what the MCP cannot show, its class in the manifest (the owner, 2026-10-15): the walk stays tested via
            # MCP, the label dictionary is not implemented there yet (v3 and v4: not_exposed on every MCP execution)
            mcp_not_implemented_check(out[-1], "label_dict_contains", LABEL_DICT_NOT_VIA_MCP)
            if key == "eco" and kind == "column":
                warm.append(region_warm_up(ctx, rec, tt, radius, radius, rec.column, sid,
                                           "the rows of a seed 868 columns carry"))
    rec = ctx.rec(LOW_CARRIER[4])
    t = ctx.flank_for(rec, ctx.window(rng, rec, 100, margin=400), 250)
    req = traverse_request([_tseed(t["window"], [rec.accession])], {"bounds": {"max_extension_bp": 200},
                                                                   "labels": {"max_labels_per_node": 1}})
    out.append(scenario(nid(), "T", "max_labels_per_node 1: the cap echoed, a cut list stated, the walk unchanged",
                        "POST /traverse", req, t,
                        [OK, INV, ex("labels_per_node", cap=1), ex("walk_in_record", label_scope="record")],
                        tags=["labels"]))
    nid = _ids("T-exhaustive")
    for i, key in enumerate(["buc", "mja", "lam"]):
        rec = ctx.rec(key)
        t = ctx.flank_for(rec, ctx.window(rng, rec, 60, margin=400), 150)
        req = traverse_request([_tseed(t["window"], [rec.accession])], {"bounds": {"max_extension_bp": 100}}, preset="trie")
        out.append(scenario(nid(), "T", "trie (exhaustive) of %s, radius 100: a leaf follows the record" % key,
                            "POST /traverse", req, t,
                            [OK, INV, ex("walk_in_record", label_scope="record")], tags=["exhaustive"]
                           ))
    # 250 bp (500 before): the trace walk read the header's coordinates at 34 bp/s cold under the load (T-coords-0003,
    # 475 of 500 bp in 15 s; T-coords-0001 took 14.6 s), within the trace preset's 15,000 ms
    nid = _ids("T-coords")
    for i, key in enumerate(["eco", "bsu", "lam"]):
        rec = ctx.rec(key)
        t = ctx.flank_for(rec, ctx.window(rng, rec, 100, margin=700), 300)
        # the trace preset (on_reconverge keep: trace support with merging is refused, SPEC §13 item 9, G15)
        req = traverse_request([_tseed(t["window"], [rec.accession])], {"bounds": {"max_extension_bp": 250}},
                               preset="trace", extra={"output": {"coordinates": True}})
        out.append(scenario(nid(), "T", "trace + coordinates on %s: positions equal the record's" % key, "POST /traverse",
                            req, t, [OK, INV, ex("trace_coordinates"), ex("walk_in_record", label_scope="record")],
                            tags=["coordinates", "trace"]))
    nid = _ids("T-budget")
    # each knob's stated stop (SPEC-labeled-traversal-core §5, §7.0): a request budget (time, memory, work) states a
    # resource_stop of its resource and may fail the seed when it cannot hold the seed itself; max_steps is a cap
    # (no resource_stop); the walk so far follows the record
    # The work budget is a deterministic cap (work is modelled): it must lie above the seed phase's cost, which only an
    # attempt's usage states (usage.per_seed[].work_seed). Its seed is therefore T-attempts-0002's eco window, whose
    # seed phase took 148,723 units on both direct arms of v4 (refseq97-6f945db4), the walk about 2,140 units per bp
    # and arm (876,572 for 200 bp on both arms): 300,000 units hold the seed phase twice over and stop the walk at a
    # few tens of bp (triage-v4 G18: a lambda seed's validation took at least 435,722 units, so 300,000 failed it on
    # every arm); the window drawn here for that slot keeps the area's stream as it was.
    SEED_PHASE_WORK = {"work_seed": 148723, "walk_units_per_bp_and_arm": 2140,
                       "source": "T-attempts-0002 usage.per_seed[0].work_seed and the arms' counters.work_units, v4 "
                                 "direct-local and direct-remote (refseq97-6f945db4)"}
    budgets = [("time_budget_ms", 1, "time", "time_budget", ["partial"]),
               ("max_memory_mb", 1, "memory", "resource_limit", ["partial", "failed"]),
               ("max_work_units", 300000, "work", "resource_limit", ["partial"]),
               ("max_steps", 5, None, "max_steps", ["partial"])]
    work_slot: Dict[str, Any] = {}
    for i, (knob, val, resource, reason, walks) in enumerate(budgets):
        rec = recs[(i * 11 + 5) % len(recs)]
        t = ctx.flank_for(rec, ctx.window(rng, rec, 100, margin=1300), 1050)
        if knob == "max_work_units":
            work_slot = {"at": len(out), "sid": nid(), "val": val, "resource": resource, "reason": reason,
                         "walks": walks}
            continue
        req = traverse_request([_tseed(t["window"], [rec.accession])], {"bounds": {"max_extension_bp": 1000, knob: val}})
        stop = {"knob": "bounds." + knob, "walks": walks, "cap_reason": reason}
        if resource:
            stop["resource"] = resource
        out.append(scenario(nid(), "T", "bounds.%s %s: its stated stop (%s), the walk so far in the record" % (
                                knob, val, resource or reason), "POST /traverse", req, t,
                            [OK, INV, ex("traverse_stop", **stop),
                             ex("walk_in_record", label_scope="record", prefix_only=True, if_walked=True)],
                            tags=["budget"]))
    nid = _ids("T-attempts")
    t = ctx.flank_for(eco, ctx.window(rng, eco, 100, margin=400), 250)
    req = traverse_request([_tseed(t["window"], [eco.accession])], {"bounds": {"max_extension_bp": 200}})
    req["not_after_ms"] = 0
    out.append(scenario(nid(), "T", "not_after_ms in the past: 409 expired, nothing run", "POST /traverse", req, t,
                        [ex("status", http=409), ex("answer_fields", equals={"state": "expired"}, any_status=True)],
                        tags=["attempts", "refusal"]))
    req = traverse_request([_tseed(t["window"], [eco.accession])], {"bounds": {"max_extension_bp": 200}})
    # an id runs once per server process while it is retained (retention_s 3600); the runner appends "-<run id>-<arm>"
    # (run_unique), so each arm of a run runs its own attempt once: T-attempts-0002 runs it (a 200 stating its usage
    # and bound), T-attempts-0003, the same request right after it, meets the retained duplicate (a 409 stating the
    # attempt). v4 asked both of one scenario's two repeats, so the 200-only checks read inconclusive on the 409
    # (triage-v4 V04) and the deploy gate needed a declared exclusion; each scenario now decides all of its checks.
    req["attempt_id"] = "suite-s1-%s" % t["window"][:16]
    # the attempt's bound (SPEC-labeled-traversal-core §6.8: seeds, the effective per-seed budget and the server's
    # allowance, capabilities attempts.allowance_ms 10,000; walk_until_ms below seeds x T + allowance)
    out.append(scenario(nid(), "T", "an attempt: the response states its usage and bound",
                        "POST /traverse", req, t,
                        [OK, INV, ex("answer_fields", present=["usage"]),
                         ex("answer_fields", equals={"usage.bound.seeds": 1, "usage.bound.time_budget_ms": 15000,
                                                     "usage.bound.allowance_ms": 10000},
                            present=["usage.bound.walk_until_ms"])],
                        tags=["attempts", "retention-dependent", "code-derived"], run_unique=["attempt_id"]))
    # the duplicate's 409: its attempt's state and usage, no response-level usage (SPEC §10.3; A-10, G17)
    out.append(scenario(nid(), "T", "the attempt id of T-attempts-0002 sent again while retained: the 409 duplicate "
                        "states its attempt", "POST /traverse", copy.deepcopy(req), t,
                        [ex("status", http=409),
                         ex("answer_fields", any_status=True, present=["attempt.state", "attempt.usage"],
                            absent=["usage"])],
                        tags=["attempts", "retention-dependent", "code-derived"], run_unique=["attempt_id"]))
    # the work budget's slot of T-budget (above): this attempt's seed, whose seed phase v4 measured
    tw = dict(ctx.flank_for(eco, t, 1050), seed_phase=SEED_PHASE_WORK)
    w = work_slot
    req = traverse_request([_tseed(tw["window"], [eco.accession])],
                           {"bounds": {"max_extension_bp": 1000, "max_work_units": w["val"]}})
    out.insert(w["at"], scenario(w["sid"], "T", "bounds.max_work_units %s above the seed phase's %d units: its stated "
                                 "stop (work), the walk so far in the record" % (w["val"], SEED_PHASE_WORK["work_seed"]),
                                 "POST /traverse", req, tw,
                                 [OK, INV, ex("traverse_stop", knob="bounds.max_work_units", walks=w["walks"],
                                              cap_reason=w["reason"], resource=w["resource"]),
                                  ex("walk_in_record", label_scope="record", prefix_only=True, if_walked=True)],
                                 tags=["budget"]))
    # T-edge
    nid = _ids("T-edge")
    t = ctx.window(rng, eco, 100, margin=400)
    base_over = {"bounds": {"max_extension_bp": 100}}
    s30 = t["window"][:30]
    out.append(scenario(nid(), "T", "seed shorter than k: refused", "POST /traverse",
                        traverse_request([_tseed(s30, [eco.accession])], base_over), t,
                        [err400("seed '%s': Seed shorter than k" % s30)], tags=["edge", "refusal", "smoke"]))
    sn = t["window"][:50] + "N" + t["window"][51:]
    out.append(scenario(nid(), "T", "seed with N: refused naming the position", "POST /traverse",
                        traverse_request([_tseed(sn, [eco.accession])], base_over), t,
                        [err400("seed '%s': Invalid character 'N' at position 50 of the seed" % sn[:32])],
                        tags=["edge", "refusal"]))
    q = st.random_dna(ctx.rng("T", "absent"), 100)
    out.append(scenario(nid(), "T", "seed absent from the graph: the request's 400 on a single-graph server", "POST /traverse",
                        traverse_request([_tseed(q, [eco.column])], base_over), {"random": True},
                        [err400(prefix="seed '%s': Seed is not fully present in the graph; graph runs: " % q[:32])],
                        tags=["edge", "negative", "refusal"]))
    for f, txt in (("graph", "Bad request: this server hosts a single graph; remove the 'graph' / 'graph_path' field"),
                   ("graphs", "Bad request: this server hosts a single graph; remove the 'graphs' field")):
        req = traverse_request([_tseed(t["window"], [eco.accession])], base_over)
        req[f] = ["x"] if f == "graphs" else "x"
        out.append(scenario(nid(), "T", "%s on a single-graph server: refused" % f, "POST /traverse", req, t,
                            [err400(txt)], tags=["edge", "refusal"]))
    out.append(scenario(nid(), "T", "a body that is not an object: refused", "POST /traverse", [t["window"]], None,
                        [err400(contains=["expected an object"])], tags=["edge", "refusal", "protocol"]))
    out.append(scenario(nid(), "T", "max_seeds + 1 (5 on staging): refused", "POST /traverse",
                        traverse_request([_tseed(ctx.window(rng, eco, 60)["window"], [eco.accession]) for _ in range(5)], base_over),
                        None, [err400("request.seeds: at most 4 seeds per request")], tags=["edge", "refusal"]))
    req = traverse_request([_tseed(t["window"], [eco.accession]), _tseed(t["window"], [eco.accession])], base_over)
    out.append(scenario(nid(), "T", "duplicate seeds: traversed again, marked duplicate", "POST /traverse", req, t,
                        [OK, INV, ex("answer_fields", equals={"results.1.duplicate": True})], tags=["edge"]))
    # the service merges identical seeds before sending, so the duplicate's result cannot exist in the MCP answer:
    # service-side by design (the owner, 2026-10-15). The table's case traverse.duplicate_seeds (mcp.by_design, from
    # the request) explains both readings the merge causes, the invariants' "1 results for 2 seeds (1 distinct)" and
    # answer_fields' "field not in the answer: results.1.duplicate", so the MCP arm runs the scenario, reads
    # service_by_design and is verified against direct-local; any other failure still fails
    req = traverse_request([_tseed(t["window"], [eco.accession])], base_over, extra={"bogus": 1})
    out.append(scenario(nid(), "T", "an unknown strategy key: refused naming it", "POST /traverse", req, t,
                        [err400("strategy: unknown field 'bogus'")], tags=["edge", "refusal"]))
    req = traverse_request([_tseed(t["window"], [eco.accession])], {"bounds": {"max_extension_bp": 100, "bogus": 1}})
    out.append(scenario(nid(), "T", "an unknown bounds key: refused naming it", "POST /traverse", req, t,
                        [err400("strategy.bounds: unknown field 'bogus'")], tags=["edge", "refusal"]))
    req = traverse_request([_tseed(t["window"], [eco.accession])], base_over, extra={"branching": {"max_label_branches": 65}})
    out.append(scenario(nid(), "T", "max_label_branches above its maximum: refused with the range", "POST /traverse", req, t,
                        [err400("strategy.branching.max_label_branches: out of range [0, 64]")], tags=["edge", "refusal"]))
    req = traverse_request([_tseed(t["window"], [eco.accession])], {"bounds": {"max_extension_bp": 100, "time_budget_ms": 70000}})
    out.append(scenario(nid(), "T", "time_budget_ms above the server's max_time_ms: lowered and stated", "POST /traverse", req, t,
                        [OK, INV, ex("clamped", where="strategy", field="bounds.time_budget_ms", requested=70000, effective=60000)],
                        tags=["edge", "clamp", "code-derived"]))
    req = traverse_request([{"sequence": t["window"], "labels": []}], base_over)
    out.append(scenario(nid(), "T", "an empty seed label list: refused", "POST /traverse", req, t,
                        [err400("request.seeds[0].labels: expected a non-empty array of labels; omit the field to derive the permitted set from the seed")],
                        tags=["edge", "refusal"]))
    req = traverse_request([_tseed(t["window"], ["NOPE_0000.1"])], base_over)
    out.append(scenario(nid(), "T", "an unknown seed label: refused naming it", "POST /traverse", req, t,
                        [err400(contains=["Label not found in the annotation: 'NOPE_0000.1'"])], tags=["edge", "refusal"]))
    req = traverse_request([_tseed(t["window"], [eco.accession])], base_over,
                           extra={"branching": {"tip_window_bp": 5}})
    out.append(scenario(nid(), "T", "tip_window_bp (not implemented): refused", "POST /traverse", req, t,
                        [err400("strategy.branching.tip_window_bp: not implemented (tip and bubble windows are not supported yet); set it to 0 or omit it")],
                        tags=["edge", "refusal"]))
    req = {"seeds": [], "strategy": traverse_request([], base_over)["strategy"]}
    out.append(scenario(nid(), "T", "no seed: refused", "POST /traverse", req, None,
                        [err400("request.seeds: expected a non-empty array")], tags=["edge", "refusal"]))
    out += _t_coverage(ctx)
    out += _t_gate_coverage(ctx, out)
    _annotate_walks(ctx, out)
    # the warm-ups, after the area's own scenarios (their ids and draws unchanged; a minimum warm-up reads the walk's
    # exact minimum, known only now); the build places each right before the walk it warms
    by_id_ = dict((x["id"], x) for x in out)
    out += warm + [w for _, sid, label, why in warm_specs for w in minimum_warm_up(ctx, by_id_[sid], label, why)]
    return out


def _t_coverage(ctx: Ctx) -> List[dict]:
    """The /traverse cases the verification found missing (M5/M8): several distinct seeds in one request, the
    graphlet detail on the direct arm, a per-arm output cap, the min_live_labels quorum, and a seed on the reverse
    strand (on this basic graph the reverse complement of a record's unique window is not in the graph). Own random
    streams, after the other T scenarios."""
    out: List[dict] = []
    rng = ctx.rng("T", "coverage")
    nid = _ids("T-seeds")
    truths = []
    for key in ("bsu", "hvo", "sso"):
        rec = ctx.rec(key)
        truths.append(ctx.flank_for(rec, ctx.window(rng, rec, 60, margin=600), 300))
    req = traverse_request([_tseed(t["window"], [t["record"]]) for t in truths],
                           {"direction": "right", "bounds": {"max_extension_bp": 200}})
    out.append(scenario(nid(), "T", "three distinct seeds of three records in one request: each walk follows its record",
                        "POST /traverse", req, truths[0],
                        [OK, INV] + [ex("walk_in_record", seed=j, label_scope="record") for j in range(3)],
                        tags=["walk", "seeds"], truths=truths))
    nid = _ids("T-graphlet")
    for key in ("mtb", "mja"):
        rec = ctx.rec(key)
        t = ctx.flank_for(rec, ctx.window(rng, rec, 60, margin=800), 550)
        req = traverse_request([_tseed(t["window"], [rec.accession])], {"bounds": {"max_extension_bp": 500}},
                               detail="graphlet")
        out.append(scenario(nid(), "T", "output.detail graphlet on %s: the MGT body's walks follow the record" % key,
                            "POST /traverse", req, t,
                            [OK, INV, ex("traverse_seed_outcome", walks=["complete", "partial"]),
                             ex("walk_in_record", label_scope="record")], tags=["walk", "graphlet"]))
    nid = _ids("T-caps")
    rec = ctx.rec("sau")
    t = ctx.flank_for(rec, ctx.window(rng, rec, 60, margin=1300), 1050)
    req = traverse_request([_tseed(t["window"], [rec.accession])],
                           {"direction": "right", "bounds": {"max_extension_bp": 1000, "max_output_bp": 50}})
    out.append(scenario(nid(), "T", "bounds.max_output_bp 50 (per arm): the arm truncated by it, its walk in the record",
                        "POST /traverse", req, t,
                        [OK, INV, ex("traverse_stop", knob="bounds.max_output_bp", cap_reason="max_output_bp",
                                     walks=["partial"]),
                         ex("walk_in_record", label_scope="record", prefix_only=True)], tags=["walk", "caps"]))
    nid = _ids("T-quorum")
    rec = ctx.rec("mpn")
    t = ctx.flank_for(rec, ctx.window(rng, rec, 60, margin=600), 300)
    req = traverse_request([_tseed(t["window"], [rec.accession])],
                           {"bounds": {"max_extension_bp": 200, "min_live_labels": 2}})
    out.append(scenario(nid(), "T", "bounds.min_live_labels 2 on a seed of one header label: every path ends "
                        "below_min_labels at once (SPEC §6.4)", "POST /traverse", req, t,
                        [OK, INV, ex("paths_end", reasons_in=["below_min_labels"], max_bp=1)],
                        tags=["labels", "quorum", "code-derived"]))
    # the record holds none of the seed's k-mers (Ctx.rc_negative_window), so the header supports none of them and the
    # seed is refused, by either of the two seed refusals: "Seed is not fully present in the graph" when a k-mer is in
    # no record of the graph, "No seed label supports every k-mer of the seed" when other records hold them all
    # (walker.cpp, traverse.cpp); which one is the index's (v3 and v4: the first, graph runs o30)
    nid = _ids("T-strand")
    rec = ctx.rec("buc")
    t = ctx.rc_negative_window(rng, rec, 60, margin=600)
    q = st.revcomp(t["window"])
    t = empirical(dict(t, strand="-", reverse_kmers_in_record=0,
                       seed_refusals=["seed '%s': Seed is not fully present in the graph; graph runs: " % q[:32],
                                      "seed '%s': No seed label supports every k-mer of the seed" % q[:32]]),
                  "which seed refusal: the reverse k-mers in no record of the graph (v3 and v4: graph runs o30)")
    out.append(scenario(nid(), "T", "the reverse complement of a buc window none of whose reverse k-mers buc holds, as "
                        "the seed under buc's header: a seed refusal",
                        "POST /traverse", traverse_request([_tseed(q, [rec.accession])], {"bounds": {"max_extension_bp": 100}}),
                        t, [err400(prefix="seed '%s': " % q[:32])],
                        tags=["edge", "negative", "reverse-complement", "refusal", "index-empirical"]))
    return out


def _branching_seed(ctx: Ctx, key: str, radius: int, tries: int = 4000) -> dict:
    """A 31-mer of |key|'s whole record with 3-8 copies on its forward strand and none reverse-complemented, whose
    header-constrained trie (manifest_lib.trie_heads: the record's own forward k-mers) holds 3 or more heads within
    radius / 2 bp on some arm: the arm a max_live_paths of 2, or a max_paths of 1, must cut before the radius. One k-mer:
    a seed phase of one row, cheap cold. Own random stream; the truth states the copies (with their flanks, for
    walk_in_record any_copy) and, per arm, the first depth with 3 heads."""
    rec = ctx.rec(key)
    seq = rec.seq
    if rec.region:
        raise ValueError("%s: a branching seed needs the whole record" % key)
    index = {side: trie_index(seq, side) for side in ("left", "right")}
    r = ctx.rng("T", "caps", key)
    for _ in range(tries):
        s = r.randrange(radius + 100, len(seq) - radius - 100)
        w = seq[s:s + K]
        if not st.is_acgt(w) or st.revcomp(w) in seq:
            continue
        copies = st.find_all(seq, w)
        if not 3 <= len(copies) <= 8:
            continue
        heads = dict((side, trie_heads(seq, s, K, side, radius, index=index[side])) for side in ("left", "right"))
        first = dict((side, next((d + 1 for d, h in enumerate(hs) if h >= 3), None)) for side, hs in heads.items())
        if not any(f is not None and f <= radius // 2 for f in first.values()):
            continue
        t = ctx.flank_for(rec, ctx.truth_at(rec, s, K), radius + 50)
        t["copies"] = [dict((f, ctx.flank_for(rec, ctx.truth_at(rec, c, K), radius + 50)[f])
                            for f in ("start", "flank", "flank_start")) for c in copies if c != s]
        t["trie"] = dict((side, {"first_depth_with_3_heads": first[side], "heads_by_depth": heads[side][:40]})
                         for side in ("left", "right"))
        return t
    raise ValueError("no branching 31-mer found in %s" % key)


def _t_gate_coverage(ctx: Ctx, built: List[dict]) -> List[dict]:
    """The traverse requests the MCP can send that no v4 scenario forced (VERIFY-D of the outside review, the gaps
    of the resolve / traverse deploy gate; P5-v4-content-06): the trie preset's per-arm caps max_live_paths and
    max_paths forced, a cap above the service's maximum (the by-design refusal of the table's
    traverse.bounds_knob_above_service_maximum against the host's answer), the permit side of min_live_labels, a
    seed of two labels one of which the host drops (seed_unsupported) with a client seed_id. Own ids (no area's
    numbering moves), own random stream or the truths of built scenarios (copied); tagged gate-coverage (gen_D leaves
    them out, so the determinism twins stay as built). The stop of a running job (traverse_stop) is the async
    session's area S."""
    out: List[dict] = []
    mge = ctx.rec("mge")
    R = 300
    t = _branching_seed(ctx, "mge", R)
    for n, (knob, val) in ((2, ("max_live_paths", 2)), (3, ("max_paths", 1))):
        req = traverse_request([_tseed(t["window"], [mge.accession])], {"bounds": {"max_extension_bp": R, knob: val}},
                               preset="trie")
        out.append(scenario("T-caps-%04d" % n, "T", "trie of a %d-copy mge 31-mer under bounds.%s %d: an arm cut by "
                            "the cap before the radius (3 heads by %s bp), its walk so far in the record" % (
                                len(t["copies"]) + 1, knob, val,
                                "/".join(str(t["trie"][x]["first_depth_with_3_heads"]) for x in ("left", "right"))),
                            "POST /traverse", req, dict(t),
                            [OK, INV, ex("traverse_stop", knob="bounds." + knob, cap_reason=knob, walks=["partial"]),
                             ex("walk_in_record", label_scope="record", any_copy=True, prefix_only=True)],
                            tags=["exhaustive", "caps", "gate-coverage"]))
    ex1 = _by_id(built, "T-exhaustive-0001")
    tb = {k: v for k, v in ex1["truth"].items() if k != "walk"}
    req = traverse_request([_tseed(tb["window"], [tb["record"]])], {"bounds": {"max_extension_bp": 100,
                                                                              "max_live_paths": 1001}}, preset="trie")
    out.append(scenario("T-caps-0004", "T", "trie with bounds.max_live_paths 1,001, above the service's maximum of "
                        "1,000: the host walks it (its own limit is none), the MCP refuses it by design",
                        "POST /traverse", req, tb, [OK, INV, ex("walk_in_record", label_scope="record")],
                        tags=["exhaustive", "caps", "gate-coverage"]))
    q1 = _by_id(built, "T-quorum-0001")
    mpn = ctx.rec(q1["truth"]["genome"])
    tq = dict(q1["truth"])
    req = traverse_request([_tseed(tq["window"], [mpn.accession, mpn.column])],
                           {"bounds": {"max_extension_bp": 200, "min_live_labels": 2}})
    out.append(scenario("T-quorum-0002", "T", "bounds.min_live_labels 2 on a seed of its header and its taxid, both "
                        "carrying every record k-mer: the walk follows the record (the permit side of T-quorum-0001)",
                        "POST /traverse", req, tq,
                        [OK, INV, ex("traverse_seed_outcome", walks=["complete", "partial"]),
                         ex("walk_in_record", label_scope="record")],
                        tags=["labels", "quorum", "gate-coverage"]))
    other = _distant(ctx, mpn, 0)
    seed_id = "suite-v5-client-seed-id"
    ts = dict(q1["truth"], dropped_label={"label": other.column, "assumes": "an unrelated taxon (another kingdom) "
                                          "lacks one of the window's k-mers, as R-multi's label_not_full does"})
    req = traverse_request([{"seed_id": seed_id, "sequence": ts["window"], "labels": [mpn.accession, other.column]}],
                           {"bounds": {"max_extension_bp": 200}})
    out.append(scenario("T-seeds-0002", "T", "a seed of two labels, %s's header and the unrelated %s (dropped: "
                        "seed_unsupported), with a client seed_id: one label left, the id echoed with its mismatch "
                        "stated, the walk in the record" % (mpn.key, other.key),
                        "POST /traverse", req, ts,
                        [OK, INV, ex("answer_fields", equals={"results.0.seed.labels_supporting_total": 1,
                                                              "results.0.seed.seed_id": seed_id,
                                                              "results.0.seed.seed_id_mismatch": True}),
                         ex("walk_in_record", label_scope="record")],
                        tags=["walk", "seeds", "multi-label", "gate-coverage"]))
    return out


def _copy_divergence(seq: str, s0: int, L: int, copies: List[int], side: str, cap: int) -> int:
    """How far the window's other copies stay identical to it on |side| (the first copy to differ makes the walk
    split there): the least common extent over the copies, at most |cap|."""
    least = cap
    for c in copies:
        if c == s0:
            continue
        n = 0
        while n < least:
            a, b = (s0 + L + n, c + L + n) if side == "right" else (s0 - 1 - n, c - 1 - n)
            if a < 0 or b < 0 or a >= len(seq) or b >= len(seq) or seq[a] != seq[b]:
                break
            n += 1
        least = min(least, n)
    return least


def _repeat_windows(ctx: Ctx, rng, n: int, ext: Optional[int] = None) -> List[Tuple[str, dict, List[int]]]:
    """Windows of 100 bp that occur at least twice (exactly) on the forward strand of a bacterial record (rRNA
    operons, IS elements): found by sampling and keeping those with 2-8 copies; with |ext|, only those whose copies
    diverge within ext - 50 bp on some side (the verification's M1: a copy identical beyond the walk's bound gives
    no branch to expect)."""
    found = []
    for key in ["eco", "bsu", "sau", "pae", "mtb", "sco"]:
        rec = ctx.rec(key)
        seq = rec.seq
        r = ctx.rng("T", "repeat", key)
        for _ in range(4000):
            s = r.randrange(2000, len(seq) - 2200)
            w = seq[s:s + 100]
            if not st.is_acgt(w) or st.revcomp(w) in seq:
                continue
            copies = st.find_all(seq, w)
            if not 2 <= len(copies) <= 8:
                continue
            t = ctx.truth_at(rec, s, 100)
            if ext is not None:
                div = {side: _copy_divergence(seq, s, 100, copies, side, ext) for side in ("left", "right")}
                if min(div.values()) > ext - 50:
                    continue
                t["copy_divergence_bp"] = div
            found.append((key, t, [rec.offset + c for c in copies]))
            break
        if len(found) >= n:
            break
    return found


def _annotate_walks(ctx: Ctx, scenarios: List[dict]) -> None:
    """truth.walk[side] = {min_bp, first, branch_by_bp} for every /traverse scenario whose walk_in_record reads
    the record's header (label_scope record) on a whole record, computed exactly from the record
    (manifest_lib.walk_expectations: one pass per strand per record). Region records (only part fetched) and
    other label scopes keep the generic minimum."""
    by_rec: Dict[str, List[Tuple[dict, dict]]] = {}
    for sc in scenarios:
        if sc["route"] != "POST /traverse" or not isinstance(sc["request"], dict):
            continue
        if not any(e["check"] == "walk_in_record" and e["args"].get("label_scope") == "record"
                   and not e["args"].get("prefix_only") for e in sc["expect"]):
            continue
        for t in sc.get("truths") or [sc["truth"]]:
            if "genome" not in t or t.get("strand", "+") != "+":
                continue
            rec = ctx.rec(t["genome"])
            if rec.region:
                continue
            by_rec.setdefault(rec.key, []).append((sc, t))
    for key, items in sorted(by_rec.items()):
        rec = ctx.rec(key)
        plans = []
        for sc, t in items:
            stg = sc["request"]["strategy"]
            d = stg.get("direction", "both")
            plans.append((t["start"] - rec.offset, t["length"], ["left", "right"] if d == "both" else [d],
                          stg["bounds"]["max_extension_bp"]))
        for (sc, t), res in zip(items, walk_expectations(rec.seq, plans)):
            t["walk"] = res


# ================================================================================================= P: /pattern
def _p(seq: str, pid: str = "w", kind: str = "dna") -> dict:
    return {"id": pid, kind: seq}


def gen_P(ctx: Ctx) -> List[dict]:
    out: List[dict] = []
    rng = ctx.rng("P")
    recs = ctx.positives
    nid = _ids("P-exact")
    i = 0
    # short patterns in suffix scope: one range, cheap; in any_offset scope a 12-mer's flank discovery ran past its
    # 60 s budget on the loaded staging host (503 deadline in the smoke run; the budget is 110,000 ms since v6,
    # PATTERN_TIME_BUDGET_MS): two such benchmark scenarios only, tagged machine-dependent
    for L, n in ((12, 4), (16, 4)):
        for j in range(n):
            rec = recs[(i * 5 + 1) % len(recs)]
            i += 1
            t = ctx.window(rng, rec, L)
            req = pattern_request([_p(t["window"])], mode="count", scope="suffix")
            out.append(scenario(nid(), "P", "%d-mer of %s in suffix scope, count: at least the record's context" % (L, rec.key),
                                "POST /pattern", req, t,
                                [OK, INV, ex("count_relation", relation="exact", min=1),
                                 ex("entry_fields", equals={"absence_scope": "suffix_only"})], tags=["positive", "count"]))
    for L, n in ((20, 5), (25, 5), (31, 5)):
        for j in range(n):
            rec = recs[(i * 5 + 1) % len(recs)]
            i += 1
            t = ctx.window(rng, rec, L)
            req = pattern_request([_p(t["window"])])
            out.append(scenario(nid(), "P", "%d-mer of %s, all_or_count: the record's contexts at their offsets" % (L, rec.key),
                                "POST /pattern", req, t,
                                [OK, INV, CC, ex("contains_source_context", strand="+", require="all"),
                                 ex("retrieval_complete", value=True)],
                                tags=["positive"] + (["smoke"] if (L == 25 and j == 0) else [])))
    nid = _ids("P-placed")
    for L in (20, 25, 31):
        for j in range(3):
            rec = recs[(i * 5 + 1) % len(recs)]
            i += 1
            t = ctx.window(rng, rec, L)
            req = pattern_request([_p(t["window"])], output={"labels": "predicate_only"},
                                  predicate={"any": [rec.column]}, predicate_strands="context")
            out.append(scenario(nid(), "P", "%d-mer of %s with any(%s), predicate_only: the record placed at the window" % (L, rec.key, rec.column),
                                "POST /pattern", req, t,
                                [OK, INV, ex("predicate_selected", selected=True, strand="+"),
                                 ex("source_placement", strand="+")],
                                tags=["positive", "placement", "predicate"] + (["smoke"] if (L == 31 and j == 0) else [])))
    for j, key in enumerate(LOW_CARRIER[:4]):
        rec = ctx.rec(key)
        t = ctx.window(rng, rec, 25 if j % 2 else 31)
        req = pattern_request([_p(t["window"])], output={"labels": "all"})
        out.append(scenario(nid(), "P", "%d-mer of %s, labels all: its taxid among the labels, the record placed" % (t["length"], key),
                            "POST /pattern", req, t,
                            [OK, INV, CC, ex("source_label", strand="+"), ex("source_placement", strand="+")],
                            tags=["positive", "placement", "labels-all"]))
    # P-scope
    nid = _ids("P-scope")
    for j in range(2):
        rec = recs[(j * 9 + 4) % len(recs)]
        t = ctx.window(rng, rec, 20)
        req = pattern_request([_p(t["window"])], scope="suffix")
        out.append(scenario(nid(), "P", "suffix scope: the record's context at offset k-L", "POST /pattern", req, t,
                            [OK, INV, CC, ex("contains_source_context", strand="+", require="all"),
                             ex("entry_fields", equals={"absence_scope": "suffix_only", "scope": "suffix"})],
                            tags=["scope"]))
    for j, key in enumerate(["lam", "t4"]):
        rec = ctx.rec(key)
        t = ctx.truth_at(rec, 0, 20)
        for scope in ("any_offset", "suffix"):
            req = pattern_request([_p(t["window"])], scope=scope)
            exp = [OK, INV, CC, ex("entry_fields", equals={"absence_scope": "any_offset" if scope == "any_offset" else "suffix_only"})]
            if scope == "any_offset":
                exp.append(ex("contains_source_context", strand="+", require="all"))
            out.append(scenario(nid(), "P", "the first 20 bases of %s in %s scope" % (key, scope), "POST /pattern",
                                req, t, exp, tags=["scope", "record-start"]))
    for j in range(2):
        rec = recs[(j * 9 + 6) % len(recs)]
        t = ctx.window(rng, rec, 25)
        req = pattern_request([_p(t["window"])], strands="forward")
        out.append(scenario(nid(), "P", "strands forward: the + contexts only", "POST /pattern", req, t,
                            [OK, INV, CC, ex("contains_source_context", strand="+", require="all"),
                             ex("entry_fields", equals={"strands": ["+"]})], tags=["strands"]))
        t2 = ctx.window(rng, rec, 25)
        t2["strand"] = "-"
        req = pattern_request([_p(st.revcomp(t2["window"]))], strands="reverse")
        out.append(scenario(nid(), "P", "strands reverse with the reverse complement: the record as a - context", "POST /pattern",
                            req, t2, [OK, INV, CC, ex("contains_source_context", strand="-", require="all"),
                                      ex("entry_fields", equals={"strands": ["-"]})], tags=["strands"]))
    # in suffix scope (one offset): the smoke validation of lane S1b saw a 16 bp palindrome's any_offset list
    # rightly withheld above 10,000 contexts after 32 s (palindromic 16-mers are enriched over chance)
    for j, (key, t) in enumerate(_palindromes(ctx, n=2)):
        t["strand"] = "="
        req = pattern_request([_p(t["window"])], scope="suffix")
        out.append(scenario(nid(), "P", "a %d bp palindrome of %s in suffix scope: searched once, strand =" % (t["length"], key),
                            "POST /pattern", req, t,
                            [OK, INV, CC, ex("entry_fields", equals={"palindromic": True, "strands": ["="]}),
                             ex("contains_source_context", strand="=", require="all")], tags=["palindrome"]))
    # P-iupac
    nid = _ids("P-iupac")
    for j in range(6):
        rec = recs[(j * 7 + 3) % len(recs)]
        t = ctx.window(rng, rec, 25)
        r2 = ctx.rng("P", "iupac", j)
        w = list(t["window"])
        npos = 2 + j % 3
        positions = sorted(r2.sample(range(25), npos))
        for p in positions:
            w[p] = st.iupac_code_holding(w[p], 2 + (p % 2), r2)
        pat = "".join(w)
        t["iupac_positions"] = positions
        req = pattern_request([_p(pat, kind="iupac")])
        out.append(scenario(nid(), "P", "25-mer with %d IUPAC codes holding the record's bases" % npos, "POST /pattern",
                            req, t, [OK, INV, CC, ex("contains_source_context", strand="+", require="all")],
                            tags=["iupac"]))
    rec = ctx.rec("bsu")
    t = ctx.window(rng, rec, 25)
    req = pattern_request([_p("NNNN" + t["window"][4:], kind="iupac")])
    out.append(scenario(nid(), "P", "a leading N run (skipped): the record still found", "POST /pattern", req, t,
                        [OK, INV, CC, ex("contains_source_context", strand="+", require="all")], tags=["iupac", "n-run"]))
    t = ctx.window(rng, rec, 25)
    req = pattern_request([_p(t["window"][:10] + "NNN" + t["window"][13:], kind="iupac")])
    out.append(scenario(nid(), "P", "an internal N run: the record still found", "POST /pattern", req, t,
                        [OK, INV, CC, ex("contains_source_context", strand="+", require="all")], tags=["iupac", "n-run"]))
    # P-modes
    nid = _ids("P-modes")
    rec = ctx.rec("sau")
    t = ctx.window(rng, rec, 25)
    out.append(scenario(nid(), "P", "count mode: no list, the count at least the record's 7 contexts", "POST /pattern",
                        pattern_request([_p(t["window"])], mode="count"), t,
                        [OK, INV, ex("count_relation", relation="exact", min=7), ex("retrieval_complete", value=False)],
                        tags=["modes"]))
    out.append(scenario(nid(), "P", "partial with max_contexts 2: cut max_contexts", "POST /pattern",
                        pattern_request([_p(t["window"])], mode="partial", max_contexts=2), t,
                        [OK, INV, CC, ex("cut", reason="max_contexts"), ex("entry_fields", equals={"returned": 2})],
                        tags=["modes"]))
    out.append(scenario(nid(), "P", "all_or_count with max_contexts 3: withheld count_above_threshold", "POST /pattern",
                        pattern_request([_p(t["window"])], max_contexts=3), t,
                        [OK, INV, CC, ex("withheld", reason="count_above_threshold"),
                         ex("count_relation", relation="exact", min=7)], tags=["modes"]))
    out.append(scenario(nid(), "P", "max_contexts above the cap: lowered and stated", "POST /pattern",
                        pattern_request([_p(t["window"])], max_contexts=10001), t,
                        [OK, INV, ex("clamped", field="max_contexts", requested=10001, effective=10000),
                         ex("contains_source_context", strand="+", require="all")], tags=["modes", "clamp"]))
    out.append(scenario(nid(), "P", "stop_at_threshold with max_contexts 3: threshold_crossed or count_above_threshold", "POST /pattern",
                        pattern_request([_p(t["window"])], max_contexts=3, stop_at_threshold=True), t,
                        [OK, INV, CC, ex("withheld", reason=["threshold_crossed", "count_above_threshold"])],
                        tags=["modes"]))
    out.append(scenario(nid(), "P", "partial with max_contexts 0: nothing returned, cut max_contexts", "POST /pattern",
                        pattern_request([_p(t["window"])], mode="partial", max_contexts=0), t,
                        [OK, INV, CC, ex("cut", reason="max_contexts"), ex("entry_fields", equals={"returned": 0})],
                        tags=["modes"]))
    # P-negative
    nid = _ids("P-negative")
    for j in range(4):
        q = st.random_dna(ctx.rng("P", "neg", j), 31)
        out.append(scenario(nid(), "P", "a random 31-mer: exact 0, complete", "POST /pattern",
                            pattern_request([_p(q)]), {"random": q},
                            [OK, INV, CC, ex("absent")], tags=["negative"]))
    for j in range(2):
        rec = recs[(j * 13 + 2) % len(recs)]
        t = ctx.window(rng, rec, 40)
        q = st.shuffled(ctx.rng("P", "shuf", j), t["window"])
        out.append(scenario(nid(), "P", "a shuffled 40-mer (paths): exact 0", "POST /pattern",
                            pattern_request([_p(q)], long_search="paths"), {"shuffled_from": t},
                            [OK, INV, CC, ex("absent", count="paths")], tags=["negative"]))
    # P-floor
    nid = _ids("P-floor")
    rec = ctx.rec("mtb")
    t12 = ctx.window(rng, rec, 12)
    t = ctx.truth_at(rec, t12["start"] - rec.offset, 14)
    w = list(t["window"])
    r2 = ctx.rng("P", "floor24")
    for p in (2, 5, 8, 11):
        w[p] = st.iupac_code_holding(w[p], 2, r2)
    q = "".join(w)
    out.append(scenario(nid(), "P", "14 positions with four two-base codes (24 bits, not exact) in suffix scope: at the floor, answered",
                        "POST /pattern", pattern_request([_p(q, kind="iupac")], mode="count", scope="suffix"), t,
                        [OK, INV, ex("entry_fields", equals={"information_bits": 24.0}, absent=["error"]),
                         ex("count_relation", relation="exact", min=1)], tags=["floor"]))
    t = ctx.window(rng, rec, 11, unique=False)
    out.append(scenario(nid(), "P", "an 11-mer in any_offset scope: below the floor", "POST /pattern",
                        pattern_request([_p(t["window"])], mode="count"), t,
                        [OK, ex("pattern_slot_error", pattern=0, code="information_below_floor",
                                text="pattern: 22.0 information bits, below the floor of 24.0 for scope any_offset")],
                        tags=["floor", "smoke"]))
    out.append(scenario(nid(), "P", "an 11-mer in suffix scope: exempt, answered", "POST /pattern",
                        pattern_request([_p(t["window"])], mode="count", scope="suffix"), t,
                        [OK, INV, ex("entry_fields", equals={"scope": "suffix"}, absent=["error"])], tags=["floor"]))
    t = ctx.window(rng, rec, 18)
    q = "NNNNNNNN" + t["window"][8:]
    out.append(scenario(nid(), "P", "N-heavy 18-mer (20 bits): below the floor", "POST /pattern",
                        pattern_request([_p(q, kind="iupac")], mode="count"), t,
                        [OK, ex("pattern_slot_error", pattern=0, code="information_below_floor",
                                text="pattern: 20.0 information bits, below the floor of 24.0 for scope any_offset")],
                        tags=["floor"]))
    t = ctx.window(rng, rec, 20)
    q = "N" * 20 + t["window"]
    out.append(scenario(nid(), "P", "N^20 + 20 bases: its forward anchor window below the floor", "POST /pattern",
                        pattern_request([_p(q, kind="iupac")], mode="count", long_search="paths"), t,
                        [OK, ex("pattern_slot_error", pattern=0, code="information_below_floor",
                                text="pattern: 22.0 information bits in the anchor window, below the floor of 24.0 for scope long")],
                        tags=["floor", "long"]))
    # P-caps
    nid = _ids("P-caps")
    rec = ctx.rec("pae")
    t = ctx.window(rng, rec, 40)
    out.append(scenario(nid(), "P", "max_anchors 0: the extension not admitted", "POST /pattern",
                        pattern_request([_p(t["window"])], long_search="paths", max_anchors=0), t,
                        [OK, INV, CC, ex("withheld", reason="anchors_above_threshold"),
                         ex("entry_fields", equals={"counts.paths.extension": "not_admitted",
                                                    "counts.paths.relation": "unknown"})], tags=["caps"]))
    t1, t2 = ctx.window(rng, rec, 25), ctx.window(rng, rec, 25)
    out.append(scenario(nid(), "P", "max_steps 1: a discovery stop, sticky for the next pattern", "POST /pattern",
                        pattern_request([_p(t1["window"], "a"), _p(t2["window"], "b")], max_steps=1), t1,
                        [OK, INV, CC, ex("stop", pattern=0, phase="discovery", reason="max_steps"),
                         ex("stop", pattern=1, reason="max_steps"), ex("withheld", pattern=0, reason="discovery_budget"),
                         ex("count_relation", pattern=0, relation=["at_least", "unknown"])],
                        tags=["caps", "smoke"]))
    out.append(scenario(nid(), "P", "max_steps 50 in partial: cut max_steps", "POST /pattern",
                        pattern_request([_p(t1["window"])], mode="partial", max_steps=50), t1,
                        [OK, INV, CC, ex("stop", phase="discovery", reason="max_steps"), ex("cut", reason="max_steps")],
                        tags=["caps"]))
    pats = [_p(ctx.window(rng, recs[j % len(recs)], 20)["window"], "p%d" % j) for j in range(16)]
    out.append(scenario(nid(), "P", "16 patterns under a 300 ms budget: each exact within it or a stated time stop "
                        "(inexact after a discovery stop, sticky), or a 503 deadline", "POST /pattern",
                        pattern_request(pats, time_budget_ms=300), None,
                        [OK, INV, CC, ex("stops_stated", reasons=["time"])],
                        tags=["caps", "machine-dependent"]))
    # P-edge
    nid = _ids("P-edge")
    t = ctx.window(rng, ctx.rec("cje"), 25)
    out.append(scenario(nid(), "P", "an empty pattern list: refused", "POST /pattern",
                        pattern_request([]), None,
                        [err400("request.patterns: expected a list of 1 to 16 pattern objects", code="invalid_request")],
                        tags=["edge", "refusal"]))
    pats = [_p(ctx.window(rng, recs[j % len(recs)], 20)["window"], "p%d" % j) for j in range(17)]
    out.append(scenario(nid(), "P", "17 patterns: refused, never cut", "POST /pattern", pattern_request(pats, mode="count"), None,
                        [err400("request.patterns: expected a list of 1 to 16 pattern objects (got 17; the server's --pattern-max-patterns)",
                                code="invalid_request")], tags=["edge", "refusal"]))
    w = t["window"]
    bad = [{"dna": w[:3] + "U" + w[4:], "id": "u"}, {"dna": w[:3] + "R" + w[4:], "id": "r"},
           {"iupac": w[:4] + "X" + w[5:], "id": "x"}, {"dna": "", "id": "empty"}, {"dna": w, "id": "ok"}]
    out.append(scenario(nid(), "P", "bad alphabet slots beside an answered pattern", "POST /pattern",
                        pattern_request(bad, mode="count"), t,
                        [OK, INV,
                         ex("pattern_slot_error", pattern=0, code="bad_alphabet", text="pattern: character 'U' at position 3 is not in the dna alphabet (A, C, G, T)"),
                         ex("pattern_slot_error", pattern=1, code="bad_alphabet", text="pattern: character 'R' at position 3 is not in the dna alphabet (A, C, G, T)"),
                         ex("pattern_slot_error", pattern=2, code="bad_alphabet", text="pattern: character 'X' at position 4 is not in the iupac alphabet (A C G T R Y S W K M B D H V N)"),
                         ex("pattern_slot_error", pattern=3, code="bad_alphabet", text="pattern: empty"),
                         ex("count_relation", pattern=4, relation="exact", min=7)],
                        tags=["edge", "slots"]))
    out.append(scenario(nid(), "P", "lower case: upper-cased, the record found", "POST /pattern",
                        pattern_request([_p(w.lower())]), t,
                        [OK, INV, ex("entry_fields", equals={"pattern": w}),
                         ex("contains_source_context", strand="+", require="all")], tags=["edge"]))
    out.append(scenario(nid(), "P", "an unknown field: refused", "POST /pattern", dict(pattern_request([_p(w)]), bogus=1), t,
                        [err400("request: unknown field 'bogus'", code="invalid_request")], tags=["edge", "refusal"]))
    out.append(scenario(nid(), "P", "duplicate ids: both answered, ids echoed", "POST /pattern",
                        pattern_request([_p(w, "same"), _p(w, "same")], mode="count"), t,
                        [OK, INV, ex("count_relation", pattern=1, relation="exact", min=7)], tags=["edge"]))
    out.append(scenario(nid(), "P", "budget_split: a later increment", "POST /pattern",
                        dict(pattern_request([_p(w)]), budget_split="even"), t,
                        [err400("request.budget_split: 'budget_split' is not served by this build", code="later_increment")],
                        tags=["edge", "refusal"]))
    out.append(scenario(nid(), "P", "graphs on a single-graph server: refused", "POST /pattern",
                        dict(pattern_request([_p(w)]), graphs=["x"]), t,
                        [err400("request.graphs: this server hosts a single graph; remove the 'graphs' field", code="invalid_request")],
                        tags=["edge", "refusal"]))
    out.append(scenario(nid(), "P", "in_ram on a single-graph server: accepted, nothing loaded", "POST /pattern",
                        dict(pattern_request([_p(w)]), in_ram=True), t,
                        [OK, INV, ex("answer_fields", equals={"timing.load_ms": 0}),
                         ex("contains_source_context", strand="+", require="all")], tags=["edge"]))
    out.append(scenario(nid(), "P", "output.occurrences true without labels: refused", "POST /pattern",
                        pattern_request([_p(w)], output={"labels": "none", "occurrences": True}), t,
                        [err400("request.output.occurrences: true needs output.labels \"all\" (occurrences are placed per label)",
                                code="invalid_request")], tags=["edge", "refusal"]))
    out.append(scenario(nid(), "P", "a pattern naming two kinds: refused", "POST /pattern",
                        pattern_request([{"id": "two", "dna": w, "iupac": w}]), t,
                        [err400("request.patterns[0]: expected exactly one of 'dna', 'iupac'", code="invalid_request")],
                        tags=["edge", "refusal"]))
    out.append(scenario(nid(), "P", "a body that is not JSON: refused", "POST /pattern", None, None,
                        [err400(code="invalid_request")], tags=["edge", "refusal", "protocol"],
                        body_raw='{"patterns": [{"dna": "%s"}],}' % w))
    # an IUPAC code at the last position, and partial with max_contexts above the count (a complete list, no cut)
    nid = _ids("P-more")
    rm = ctx.rng("P", "more")
    rec = ctx.rec("mja")
    t = ctx.window(rm, rec, 25)
    w = list(t["window"])
    w[24] = st.iupac_code_holding(w[24], 2, rm)
    w[0] = st.iupac_code_holding(w[0], 3, rm)
    t["iupac_positions"] = [0, 24]
    out.append(scenario(nid(), "P", "IUPAC codes at the first and the last position: the record's contexts found",
                        "POST /pattern", pattern_request([_p("".join(w), kind="iupac")]), t,
                        [OK, INV, CC, ex("contains_source_context", strand="+", require="all")], tags=["iupac"]))
    rec = ctx.rec("hvo")
    t = ctx.window(rm, rec, 25)
    out.append(scenario(nid(), "P", "partial with max_contexts 10,000 above the count: the complete list, no cut",
                        "POST /pattern", pattern_request([_p(t["window"])], mode="partial"), t,
                        [OK, INV, CC, ex("contains_source_context", strand="+", require="all"),
                         ex("retrieval_complete", value=True), ex("entry_fields", equals={"cut": None})],
                        tags=["modes"]))
    # P-heavy: short patterns in any_offset scope. A 12-mer's flank discovery ran past its 60 s budget on the
    # loaded staging host (503 deadline, three times, in the smoke run; at ~60.2 s again in the v5 baseline run,
    # hence v6's 110,000 ms, PATTERN_TIME_BUDGET_MS): kept as benchmark scenarios, a 200 with its counts or the
    # stated 503 (an own random stream, after the others, so that their windows stay as built)
    nid = _ids("P-heavy")
    rh = ctx.rng("P", "heavy")
    for j, L in enumerate((12, 16)):
        rec = recs[(j * 5 + 3) % len(recs)]
        t = ctx.window(rh, rec, L)
        nsrc = len(source_contexts(t))
        out.append(scenario(nid(), "P", "%d-mer of %s in any_offset scope, count: at least its %d record contexts (heavy)" % (L, rec.key, nsrc),
                            "POST /pattern", pattern_request([_p(t["window"])], mode="count"), t,
                            [OK, INV, ex("count_relation", relation=["exact", "at_least"], min=nsrc)],
                            timing_class="slow", tags=["positive", "count", "machine-dependent", "heavy"]))
    return out


def _palindromes(ctx: Ctx, n: int) -> List[Tuple[str, dict]]:
    out = []
    # 16 bp only: a 14-bp palindrome meets ~42,000 chance contexts in any_offset scope, above the cap of 10,000
    # (the verification's S1); the scenarios search them in suffix scope (a single offset) besides
    for key in ["lam", "t4", "eco", "bsu"]:
        rec = ctx.rec(key)
        seq = rec.seq
        for L in (16,):
            r = ctx.rng("P", "pal", key, L)
            starts = list(range(100, len(seq) - 200))
            r.shuffle(starts)
            for s in starts[:200000]:
                w = seq[s:s + L]
                if w == st.revcomp(w) and st.is_acgt(w) and st.low_complexity_dust_score(w) < 1.5:
                    out.append((key, ctx.truth_at(rec, s, L)))
                    break
            if out and out[-1][0] == key:
                break
        if len(out) >= n:
            break
    return out


# ================================================================================================= L: long patterns
def gen_L(ctx: Ctx) -> List[dict]:
    out: List[dict] = []
    rng = ctx.rng("L")
    recs = ctx.positives
    nid = _ids("L-paths")
    i = 0
    for L in (32, 40, 60, 100, 200):
        for j in range(3):
            rec = recs[(i * 7 + 5) % len(recs)]
            i += 1
            t = ctx.window(rng, rec, L)
            req = pattern_request([_p(t["window"])], long_search="paths")
            out.append(scenario(nid(), "L", "%d bp of %s with long_search paths: a path spells it" % (L, rec.key),
                                "POST /pattern", req, t,
                                [OK, INV, CC, ex("path_spells", strand="+"), ex("count_relation", count="paths", relation="exact", min=1),
                                 ex("entry_fields", equals={"counts.paths.extension": "completed", "scope": "long"})],
                                tags=["paths"] + (["smoke"] if (L == 60 and j == 0) else [])))
    nid = _ids("L-verified")
    for j, key in enumerate(LOW_CARRIER[:6]):
        rec = ctx.rec(key)
        t = ctx.window(rng, rec, [40, 60, 100][j % 3])
        req = pattern_request([_p(t["window"])], long_search="paths", output={"labels": "all"},
                              require_support="record_verified")
        out.append(scenario(nid(), "L", "%d bp of %s, paths with labels, require_support record_verified" % (t["length"], key),
                            "POST /pattern", req, t,
                            [OK, INV, CC, ex("path_record_verified", strand="+")], tags=["paths", "labels-all", "record-verified"]))
    nid = _ids("L-supported")
    for j in range(8):
        rec = recs[(j * 5 + 2) % len(recs)]
        t = ctx.window(rng, rec, [40, 60, 100, 200][j % 4])
        req = pattern_request([_p(t["window"])], long_search="supported_paths", output={"labels": "predicate_only"},
                              predicate={"any": [rec.column]}, predicate_strands="context", max_paths=1024)
        out.append(scenario(nid(), "L", "%d bp of %s: supported paths selected by any(%s), the record verified" % (t["length"], rec.key, rec.column),
                            "POST /pattern", req, t,
                            [OK, INV, ex("path_record_verified", strand="+"),
                             ex("count_relation", count="supported_paths", min=1)],
                            tags=["supported-paths", "predicate"] + (["smoke"] if j == 0 else [])))
    for j in range(2):
        rec = recs[(j * 5 + 9) % len(recs)]
        t = ctx.window(rng, rec, 60)
        req = pattern_request([_p(t["window"])], long_search="supported_paths", supported_paths_level="label_intersection",
                              max_paths=1024)
        out.append(scenario(nid(), "L", "supported paths at the label_intersection level", "POST /pattern", req, t,
                            [OK, INV, ex("path_spells", strand="+"),
                             ex("entry_fields", equals={"counts.supported_paths.level": "label_intersection"})],
                            tags=["supported-paths"]))
    nid = _ids("L-anchors")
    for j in range(3):
        rec = recs[(j * 11 + 1) % len(recs)]
        t = ctx.window(rng, rec, 45)
        req = pattern_request([_p(t["window"])])
        out.append(scenario(nid(), "L", "45 bp without long_search: anchors counted, paths withheld", "POST /pattern", req, t,
                            [OK, INV, ex("withheld", reason="paths_later_increment"), ex("notes", has=["paths_later_increment"]),
                             ex("count_relation", count="anchors", relation="exact", min=1),
                             ex("entry_fields", equals={"counts.paths.relation": "unknown"})], tags=["anchors"]))
    nid = _ids("L-caps")
    rec = ctx.rec("hvo")
    t = ctx.window(rng, rec, 50)
    out.append(scenario(nid(), "L", "partial with max_paths 0: cut max_paths", "POST /pattern",
                        pattern_request([_p(t["window"])], mode="partial", long_search="paths", max_paths=0), t,
                        [OK, INV, ex("cut", reason="max_paths"), ex("entry_fields", equals={"returned": 0})], tags=["paths", "caps"]))
    out.append(scenario(nid(), "L", "all_or_count with max_paths 0: withheld count_above_threshold", "POST /pattern",
                        pattern_request([_p(t["window"])], long_search="paths", max_paths=0), t,
                        [OK, INV, ex("withheld", reason="count_above_threshold")], tags=["paths", "caps"]))
    out.append(scenario(nid(), "L", "supported paths with max_annotation_work 10: a stated stop of the search", "POST /pattern",
                        pattern_request([_p(t["window"])], long_search="supported_paths", max_annotation_work=10), t,
                        [OK, INV, ex("withheld", reason="annotation_budget"),
                         ex("stop", phase="extension", reason="max_annotation_work")], tags=["supported-paths", "caps"]))
    out.append(scenario(nid(), "L", "supported paths with max_annotation_work 10 in partial: cut", "POST /pattern",
                        pattern_request([_p(t["window"])], mode="partial", long_search="supported_paths", max_annotation_work=10), t,
                        [OK, INV, ex("cut", reason="max_annotation_work")], tags=["supported-paths", "caps"]))
    nid = _ids("L-memory")
    out.append(scenario(nid(), "L", "A^100000 on a self-loop, paths, max_memory_mb 1: stop {extension, max_memory}", "POST /pattern",
                        pattern_request([_p("A" * 100000, "polyA")], mode="partial", long_search="paths", max_memory_mb=1), None,
                        [OK, INV, ex("stop", phase=["extension", "output"], reason="max_memory"), ex("cut", reason="max_memory"),
                         ex("notes", lacks=["low_complexity_pattern"])],
                        timing_class="slow", tags=["paths", "memory", "low-complexity"]))
    # marker-gene lengths (the verification's M5): long windows of low-carrier records as paths; one on the minus
    # strand; an IUPAC long pattern; a max_anchors threshold below the anchors (the extension not admitted)
    nid = _ids("L-long")
    rl = ctx.rng("L", "long")
    for j, (key, L) in enumerate((("buc", 500), ("mja", 1000), ("cal", 2000), ("hvo", 500))):
        rec = ctx.rec(key)
        t = ctx.window(rl, rec, L, max_dust=3.0)
        out.append(scenario(nid(), "L", "%d bp of %s with long_search paths: one path spells it" % (L, key),
                            "POST /pattern", pattern_request([_p(t["window"])], long_search="paths"), t,
                            [OK, INV, CC, ex("path_spells", strand="+"),
                             ex("count_relation", count="paths", relation="exact", min=1)],
                            timing_class="slow" if L >= 1000 else "fast", tags=["paths", "long-window"]))
    rec = ctx.rec("mpn")
    t = ctx.window(rl, rec, 100)
    t["strand"] = "-"
    # the record's path listed as '-' (P5-v4-content-04: the invariants already pin an exact pattern's '-' paths to the
    # record's bases, VERIFY-C; the title's claim is now a check of its own)
    out.append(scenario(nid(), "L", "100 bp of mpn as its reverse complement, strands reverse: the record's path as '-'",
                        "POST /pattern", pattern_request([_p(st.revcomp(t["window"]))], long_search="paths",
                                                         strands="reverse"), t,
                        [OK, INV, CC, ex("count_relation", count="paths", relation="exact", min=1),
                         ex("path_spells", strand="-")],
                        tags=["paths", "strands"]))
    rec = ctx.rec("sso")
    t = ctx.window(rl, rec, 60)
    w = list(t["window"])
    for p in (10, 30, 50):
        w[p] = st.iupac_code_holding(w[p], 2, rl)
    t["iupac_positions"] = [10, 30, 50]
    # the record's own path, on '+' (finding 6A of the outside review, VERIFY-C: a count of 1 or more passed with the
    # planted '+' path dropped and only the '-' instance listed)
    out.append(scenario(nid(), "L", "a 60 bp IUPAC pattern (3 two-base codes) with paths: the record's path among them",
                        "POST /pattern", pattern_request([_p("".join(w), kind="iupac")], long_search="paths"), t,
                        [OK, INV, CC, ex("count_relation", count="paths", relation="exact", min=1),
                         ex("path_spells", strand="+")],
                        tags=["paths", "iupac"]))
    # an exact anchor window is one k-mer per orientation: A^60's anchors are A^31 (+) and T^31 (-), two, above
    # max_anchors 1
    out.append(scenario(nid(), "L", "A^60 with paths and max_anchors 1: its two anchors (A^31, T^31) above it, the "
                        "extension not admitted", "POST /pattern",
                        pattern_request([_p("A" * 60, "homo")], long_search="paths", max_anchors=1), {"low_complexity": "A^60"},
                        [OK, INV, CC, ex("withheld", reason="anchors_above_threshold"),
                         ex("count_relation", count="anchors", relation="exact", min=2, max=2),
                         ex("entry_fields", equals={"counts.paths.extension": "not_admitted"})],
                        tags=["paths", "caps", "low-complexity", "code-derived"]))
    nid = _ids("L-repeat")
    for j, q in enumerate(["A" * 60, "T" * 100]):
        out.append(scenario(nid(), "L", "homopolymer of %d with paths: low_complexity_pattern noted" % len(q), "POST /pattern",
                            pattern_request([_p(q, "homo")], long_search="paths"), None,
                            [OK, INV, CC, ex("notes", has=["low_complexity_pattern"]),
                             ex("count_relation", count="paths", relation="exact", min=1)], tags=["paths", "low-complexity"]))
    return out


# ================================================================================================= X: peptides
def _peptide_truth(ctx: Ctx, prot: dict, win: dict, dna_override: Optional[str] = None) -> Optional[dict]:
    rec = ctx.rec(prot["genome"])
    if not rec.usable_for_positives:
        return None
    blocks = win["genomic_blocks"]
    if len(blocks) != 1:
        return None
    a, b = blocks[0]                      # 1-based inclusive on the record's forward strand
    s = a - 1 - rec.offset
    if s < 0 or b - rec.offset > len(rec.seq):
        return None
    L = b - a + 1
    t = ctx.truth_at(rec, s, L)
    strand = "+" if prot["strand"] == 1 else "-"
    dna = win["dna"]
    if (strand == "+" and t["window"] != dna) or (strand == "-" and t["window"] != st.revcomp(dna)):
        return None                       # not the record's bases: never used
    t["strand"] = strand
    t["peptide"] = win["peptide"]
    t["coding_dna"] = dna
    t["protein"] = prot["key"]
    t["transl_table"] = prot["transl_table"]
    t["window_id"] = win["id"]
    return t


def _peptide_windows(prot: dict, length_aa: int) -> List[dict]:
    """The manifest_v1 windows of |length_aa| residues of a protein, its shifted windows included."""
    out = []
    for win in prot["windows"]:
        if win["length_aa"] != length_aa:
            continue
        out.append(win)
        for sh in win.get("shifts") or []:
            out.append(dict(win, id=sh["id"], peptide=sh["peptide"], dna=sh["dna"], genomic_blocks=sh["genomic_blocks"]))
    return out


def _table_positives(ctx: Ctx, prots: List[dict], length_aa: int, quota: Dict[int, int]) -> List[Tuple[dict, dict, dict]]:
    """Per genetic code, up to quota[table] windows of |length_aa| residues taken round-robin over the code's
    genomes; for a code other than 1 and 11, windows whose coding DNA reads differently in table 1 first (the
    verification's M4: else a positive proves nothing about the table). Returns (protein, window, truth)."""
    out = []
    for table, n in sorted(quota.items()):
        by_genome: Dict[str, List[Tuple[int, dict, dict, dict]]] = {}
        for prot in prots:
            if prot["transl_table"] != table:
                continue
            for win in _peptide_windows(prot, length_aa):
                t = _peptide_truth(ctx, prot, win)
                if t is None:
                    continue
                differs = st.translate(win["dna"], 1) != win["peptide"]
                rank = 0 if (table in (1, 11) or differs) else 1
                t["table1_reads_differently"] = differs
                by_genome.setdefault(prot["genome"], []).append((rank, prot, win, t))
        for g in by_genome.values():
            g.sort(key=lambda x: x[0])
        picked: List[Tuple[dict, dict, dict]] = []
        while len(picked) < n and any(by_genome.values()):
            best = min((g[0][0], k) for k, g in by_genome.items() if g)
            for k in sorted(by_genome):
                g = by_genome[k]
                if g and g[0][0] == best[0] and len(picked) < n:
                    _, prot, win, t = g.pop(0)
                    picked.append((prot, win, t))
        if len(picked) < n:
            ctx.skip("X", "%d-residue positives of table %d" % (length_aa, table),
                     "%d of the %d planned exist in the fetched records" % (len(picked), n))
        alike = sum(1 for _, _, t in picked if not t["table1_reads_differently"])
        if table not in (1, 11) and alike:
            ctx.skip("X", "%d-residue positives of table %d read differently in table 1" % (length_aa, table),
                     "%d of %d windows hold no codon the table reassigns (all its genomes' windows and shifts "
                     "searched): they test the table's sense codons only" % (alike, len(picked)))
        out += picked
    return out


def gen_X(ctx: Ctx) -> List[dict]:
    from truth import load_proteins
    out: List[dict] = []
    prots = [p for p in load_proteins() if p["genome"] in ctx.by_key]
    nid10 = _ids("X-peptide10")
    nidl = _ids("X-peptide-long")
    # every genetic code the genomes carry (1, 2, 3, 4, 5, 6, 11, 12), balanced per code
    quota10 = {1: 8, 2: 5, 3: 5, 4: 5, 5: 5, 6: 5, 11: 8, 12: 5}
    for n10, (prot, win, t) in enumerate(_table_positives(ctx, prots, 10, quota10)):
        table = prot["transl_table"]
        req = pattern_request([{"id": win["id"], "protein": win["peptide"]}], genetic_code=table)
        out.append(scenario(nid10(), "X", "decapeptide %s of %s (table %d%s): the coding DNA's contexts" % (
                                win["id"], prot["key"], table,
                                ", read differently in table 1" if t["table1_reads_differently"] else ""),
                            "POST /pattern", req, t,
                            [OK, INV, CC, ex("contains_source_context", strand=t["strand"], require="all")],
                            tags=["peptide", "table-%d" % table] + (["smoke"] if n10 == 0 else [])))
    quotal = {1: 4, 2: 3, 3: 3, 4: 3, 5: 3, 6: 3, 11: 4, 12: 3}
    picked = _table_positives(ctx, prots, 20, quotal)
    picked += [x for x in _table_positives(ctx, prots, 30, {1: 1, 11: 1})]
    for prot, win, t in picked:
        table = prot["transl_table"]
        req = pattern_request([{"id": win["id"], "protein": win["peptide"]}], genetic_code=table,
                              long_search="supported_paths", output={"labels": "predicate_only"},
                              predicate={"any": [t["taxid"]]}, predicate_strands="context", max_paths=1024)
        out.append(scenario(nidl(), "X", "%d-residue %s of %s (table %d): a supported path of its taxid, record verified" % (
                                win["length_aa"], win["id"], prot["key"], table),
                            "POST /pattern", req, t,
                            [OK, INV, ex("path_record_verified", strand=t["strand"])],
                            timing_class="slow", tags=["peptide", "supported-paths", "table-%d" % table]))
    # the table-1 arm of non-standard codes: where a reassigned codon lies in the window the source is out of reach
    nid = _ids("X-table1")
    n = 0
    built_tables: set = set()
    for prot in prots:
        if prot["transl_table"] in (1, 11):
            continue
        for win in prot["windows"]:
            if win["length_aa"] != 10:
                continue
            t = _peptide_truth(ctx, prot, win)
            if t is None:
                continue
            t1 = st.translate(t["coding_dna"], 1)
            if t1 == win["peptide"]:
                continue   # no reassigned codon: table 1 reads it alike
            req = pattern_request([{"id": win["id"], "protein": win["peptide"]}], genetic_code=1)
            out.append(scenario(nid(), "X", "decapeptide of %s (table %d) read in table 1: the record's codons are out of reach" % (prot["genome"], prot["transl_table"]),
                                "POST /pattern", req, t,
                                [OK, INV, CC, ex("source_absent", strand=t["strand"])], tags=["peptide", "wrong-code", "expected-miss"]))
            built_tables.add(prot["transl_table"])
            n += 1
            break
        if n >= 6:
            break
    # every other code the genomes carry without a scenario is stated (P5-v4-content-05: nothing drops silently)
    for table in sorted(set(p["transl_table"] for p in prots) - {1, 11} - built_tables):
        reads = [w for p in prots if p["transl_table"] == table for w in p["windows"]
                 if w["length_aa"] == 10 and st.translate(w["dna"], 1) != w["peptide"]]
        ctx.skip("X", "X-table1 table %d" % table,
                 ("the six scenarios are taken from the first proteins with such a window (tables %s); table %d's "
                  "%d windows that read differently in table 1 were not reached" % (
                      ", ".join(str(x) for x in sorted(built_tables)), table, len(reads))) if reads else
                 "no decapeptide of table %d reads differently in table 1 (no codon the table reassigns)" % table)
    # X-stop: the last 9 residues and the stop codon
    nid = _ids("X-stop")
    n = 0
    per_table: Dict[int, int] = {}
    for prot in prots:
        # a stop positive per code with stop codons that the genomes carry: 1, 2, 4, 11, 12 (one each, two for 11)
        if prot["transl_table"] not in (1, 2, 4, 11, 12) or prot["n_exons"] != 1 \
                or not ctx.rec(prot["genome"]).usable_for_positives \
                or per_table.get(prot["transl_table"], 0) >= (2 if prot["transl_table"] == 11 else 1):
            continue
        rec = ctx.rec(prot["genome"])
        cds = prot["cds_seq"]
        tail = cds[-30:]
        pep = prot["protein_seq"][-9:] + "*"
        if st.translate(tail, prot["transl_table"]) != pep:
            continue
        a, b = prot["spans"][0]
        if prot["strand"] == 1:
            s = b - 30 - rec.offset
        else:
            s = a - 1 - rec.offset
        if s < 0 or s + 30 > len(rec.seq):
            continue
        t = ctx.truth_at(rec, s, 30)
        exp_bases = tail if prot["strand"] == 1 else st.revcomp(tail)
        if t["window"] != exp_bases:
            continue
        t["strand"] = "+" if prot["strand"] == 1 else "-"
        t["peptide"] = pep
        req = pattern_request([{"id": "stop", "protein": pep}], genetic_code=prot["transl_table"])
        out.append(scenario(nid(), "X", "%s's last 9 residues and its stop (table %d): found at the stop codon" % (prot["key"], prot["transl_table"]),
                            "POST /pattern", req, t,
                            [OK, INV, CC, ex("contains_source_context", strand=t["strand"], require="all")],
                            tags=["peptide", "stop", "table-%d" % prot["transl_table"]]))
        per_table[prot["transl_table"]] = per_table.get(prot["transl_table"], 0) + 1
        n += 1
    for table in (1, 2, 4, 11, 12):
        if not per_table.get(table):
            ctx.skip("X", "X-stop table %d" % table, "no single-exon protein whose last 9 residues and stop read in "
                     "table %d lie in a fetched record" % table)
    for table in (27, 28, 31):
        req = pattern_request([{"id": "nostop", "protein": "MKV*LLAAGR"}], genetic_code=table)
        out.append(scenario(nid(), "X", "a stop in table %d (no stop codon): exact 0 without a search" % table, "POST /pattern",
                            req, None,
                            [OK, INV, ex("notes", has=["no_stop_codon"]), ex("absent"), ex("stop", null=True),
                             ex("entry_fields", equals={"work.steps": 0, "determinism": "full"})],
                            tags=["peptide", "stop"] + (["smoke"] if table == 27 else [])))
    # X-junction: expected 0 paths spelling the spliced coding DNA in its record
    nid = _ids("X-junction")
    n = 0
    for prot in prots:
        for jw in prot.get("junction_windows") or []:
            rec = ctx.rec(prot["genome"])
            if not rec.usable_for_positives or "dna" not in jw:
                continue
            seq = jw["dna"]
            t = empirical({"record": rec.accession, "taxid": rec.column,
                           "spliced_dna": seq if prot["strand"] == 1 else st.revcomp(seq),
                           "peptide": jw["peptide"], "protein": prot["key"], "window_id": jw.get("id")},
                          "no path of the whole graph spells the spliced DNA: its record holds it on neither strand, "
                          "another record (an intronless copy, a cDNA) or a recombination of the graph could")
            req = pattern_request([{"id": jw.get("id", "junction"), "protein": jw["peptide"]}], genetic_code=prot["transl_table"],
                                  mode="partial", long_search="paths", max_paths=1024)
            out.append(scenario(nid(), "X", "window across one exon junction of %s: no path spells the spliced DNA" % prot["key"],
                                "POST /pattern", req, t, [OK, INV, ex("path_absent")],
                                timing_class="slow", tags=["peptide", "junction", "expected-miss", "index-empirical"]))
            n += 1
            break
        if n >= 4:
            break
    # X-edge
    nid = _ids("X-edge")
    p0 = next(p for p in prots if p["transl_table"] == 11 and ctx.rec(p["genome"]).usable_for_positives)
    w0 = next(w for w in p0["windows"] if w["length_aa"] == 10)
    t = _peptide_truth(ctx, p0, w0)
    pep = w0["peptide"]
    # B (D|N), Z (E|Q), J (I|L) and X each in place of a member residue of a table-11 decapeptide that holds one
    for code, members in (("B", "DN"), ("Z", "EQ"), ("J", "IL"), ("X", "")):
        cand = None
        for prot in prots:
            if prot["transl_table"] != 11 or not ctx.rec(prot["genome"]).usable_for_positives:
                continue
            for win in _peptide_windows(prot, 10):
                pos = [k for k, aa in enumerate(win["peptide"]) if (aa in members if members else 3 <= k <= 6)]
                tt = _peptide_truth(ctx, prot, win) if pos else None
                if tt is not None:
                    cand = (win, tt, pos[0])
                    break
            if cand:
                break
        if cand is None:
            ctx.skip("X", "X-edge %s" % code, "no table-11 decapeptide holds a member of %s" % code)
            continue
        win, tt, k = cand
        v = win["peptide"][:k] + code + win["peptide"][k + 1:]
        tt["ambiguity"] = {"residue": k, "code": code, "replaces": win["peptide"][k]}
        req = pattern_request([{"id": code, "protein": v}], genetic_code=11)
        out.append(scenario(nid(), "X", "a decapeptide with %s for its residue %s: still found" % (code, win["peptide"][k]),
                            "POST /pattern", req, tt,
                            [OK, INV, CC, ex("contains_source_context", strand=tt["strand"], require="all")],
                            tags=["peptide", "ambiguity"]))
    req = pattern_request([{"id": "gc7", "protein": pep}], genetic_code=7)
    out.append(scenario(nid(), "X", "genetic_code 7: refused", "POST /pattern", req, t,
                        [err400("request.genetic_code: 7 is not an NCBI translation table id; the tables served are 1-6, 9-16, 21-33 (capabilities genetic_codes), 1, the standard code, the default",
                                code="genetic_code_unknown")], tags=["peptide", "refusal"]))
    w11 = next(w for w in p0["windows"] if w["length_aa"] == 20)["peptide"][:11]
    req = pattern_request([{"id": "p11", "protein": w11}], genetic_code=11)
    out.append(scenario(nid(), "X", "11 residues (33 bases > k) without long_search: anchors only", "POST /pattern", req, None,
                        [OK, INV, ex("withheld", reason="paths_later_increment"),
                         ex("entry_fields", equals={"length": 33, "residues": 11, "scope": "long"})], tags=["peptide", "long"]))
    req = pattern_request([{"id": "u", "protein": pep[:4] + "U" + pep[5:]}, {"id": "ok", "protein": pep}], genetic_code=11)
    out.append(scenario(nid(), "X", "selenocysteine U: bad_alphabet slot, the other answered", "POST /pattern", req, t,
                        [OK, INV, ex("pattern_slot_error", pattern=0, code="bad_alphabet",
                                     text="pattern: character 'U' at position 4 is not in the protein alphabet (A C D E F G H I K L M N P Q R S T V W Y, and X B Z J)"),
                         ex("contains_source_context", pattern=1, strand=t["strand"], require="all")], tags=["peptide", "slots"]))
    return out


# ================================================================================================= M: mismatches
def _mm_window(ctx: Ctx, rng, rec, L: int) -> dict:
    return ctx.window(rng, rec, L)


REFSEQ33M_KMERS = 626.2e9


def _expected_ball_contexts(L: int, m: int) -> float:
    """Contexts a random exact L-mer's ball of m substitutions meets by chance on refseq33m (both strands, every
    offset): |B| x k-mers x (k - L + 1) x 2 / 4^L."""
    return st.mismatch_ball_size_iupac("A" * L, m) * REFSEQ33M_KMERS * (K - L + 1) * 2 / 4 ** L


# A planted probe's list must be released whole (all_or_count withholds a list above max_contexts 10,000): its length
# is the shortest candidate whose chance estimate stays below MM_CONTEXTS_TARGET, a tenth of max_contexts. The
# uniform model above is blind to composition: M-planted-0001 (a 75 % AT 20-mer of M. genitalium at m 1, 834 by the
# model) met 16,264 contexts (G03, as VERIFY-MN has it). The estimate is therefore scaled by the window's likelihood
# under its record's base composition against the uniform one (prod over its bases of 4 f(b): 11.3 for that window,
# 9,418 contexts), where that factor exceeds 1: a window less likely than uniform is not credited (the run's windows
# with a factor below 1 met up to 15 times their scaled estimate). The target's tenfold margin under max_contexts is
# the safety factor: M-planted-0001 met 1.7 times its scaled estimate, and no planted window of the run estimated
# above 50 contexts met more than 3.6 times it.
MM_CONTEXTS_TARGET = 1000.0
_COMPOSITION: Dict[str, Dict[str, float]] = {}


def _composition(rec) -> Dict[str, float]:
    """The base frequencies of a record's fetched sequence (ACGT only)."""
    if rec.key not in _COMPOSITION:
        n = {b: rec.seq.count(b) for b in "ACGT"}
        tot = float(sum(n.values())) or 1.0
        _COMPOSITION[rec.key] = {b: n[b] / tot for b in "ACGT"}
    return _COMPOSITION[rec.key]


def _composition_factor(rec, window: str) -> float:
    """The window's likelihood under its record's composition against the uniform one: prod of 4 f(b)."""
    f = _composition(rec)
    x = 1.0
    for b in window:
        x *= 4.0 * f.get(b, 0.25)
    return x


def _planted_window(ctx: Ctx, rng, rec, candidates: List[int], m: int, redraw_rng) -> dict:
    """The window of the shortest candidate length whose chance estimate, the uniform ball model scaled by the
    window's composition factor where it exceeds 1, stays below MM_CONTEXTS_TARGET (a candidate whose uniform
    estimate alone does not is passed over before a window is drawn); the last candidate is taken whatever its
    estimate. The first window is drawn from |rng| as before the composition was read, a window drawn after a
    rejected one from |redraw_rng|, so that the area's other windows keep theirs. The truth states the estimate
    (chance_contexts_estimate). After the first candidate only longer ones are tried, shortest first, so the last
    one taken whatever its estimate is the longest (A-06: [31, 25, 31] tried 25 after 31)."""
    cands = [candidates[0]] + sorted(set(x for x in candidates[1:] if x > candidates[0]))
    drawn = 0
    for n, L in enumerate(cands):
        last = n == len(cands) - 1
        base = _expected_ball_contexts(L, m)
        if not last and base >= MM_CONTEXTS_TARGET:
            continue
        tr = _mm_window(ctx, redraw_rng if drawn else rng, rec, L)
        drawn += 1
        est = base * max(1.0, _composition_factor(rec, tr["window"]))
        if est < MM_CONTEXTS_TARGET or last:
            tr["chance_contexts_estimate"] = round(est, 1)
            return tr
    raise AssertionError("no candidate length")


def _subs_plus(variant: str, record_bases: str, positions: List[int]) -> List[dict]:
    return [{"position": p, "expected": variant[p], "found": record_bases[p]} for p in sorted(positions)]


def _subs_minus(pattern: str, record_bases: str) -> List[dict]:
    """P on the reverse orientation: instance = record bases; positions in P's coordinates (SPEC §21.8)."""
    L = len(pattern)
    read = st.revcomp(record_bases)
    return [{"position": i, "expected": pattern[i], "found": read[i]} for i in range(L) if read[i] != pattern[i]]


def _above_kind_cap(req: Any) -> Optional[str]:
    """Why a /pattern request's max_mismatches lies above the snapshot's cap of one of its patterns' kinds (a pattern
    of at most k bases; manifest_lib.HOST_MISMATCH_MAX), else None."""
    m = req.get("max_mismatches") if isinstance(req, dict) else None
    if not isinstance(m, int) or isinstance(m, bool) or m < 1:
        return None
    for kind, n in _pattern_kinds(req.get("patterns")):
        cap = HOST_MISMATCH_MAX.get(kind, 3)
        if n <= K and m > cap:
            return "max_mismatches %d above the snapshot's %s cap of %d" % (m, kind, cap)
    return None


def _within_kind_caps(ctx: Ctx, area: str, scenarios: List[dict]) -> List[dict]:
    """|scenarios| without those whose max_mismatches the snapshot's cap of a pattern kind does not admit (A-07:
    above it the host answers the pattern's slot mismatches_above_cap, which no positive, negative, floor or order
    expectation states), each listed in `skipped` (its id consumed): a scenario that expects the slot refusal, or a
    refusal of the request (a status that admits no 200: manifest_lib.allowed_statuses, SD-094; [200, 400] admits
    the answer), stays."""
    out = []
    for s in scenarios:
        why = _above_kind_cap(s.get("request")) if s.get("route") == "POST /pattern" else None
        refusal = 200 not in allowed_statuses(s.get("expect"))
        if why and not refusal and "mismatches_above_cap" not in json.dumps(s.get("expect")):
            ctx.skip(area, "%s (%s)" % (s["id"], s["title"]), why + ": the host answers its slot, not the expectation")
            continue
        out.append(s)
    return out


def gen_M(ctx: Ctx) -> List[dict]:
    return _within_kind_caps(ctx, "M", _gen_M(ctx))


def _gen_M(ctx: Ctx) -> List[dict]:
    out: List[dict] = []
    rng = ctx.rng("M")
    # planted windows from records with few relatives in RefSeq (the verification's M2: a window of a gene
    # conserved across many genera has homologs within 1-3 substitutions that can push the list past the cap;
    # bacteria, archaea and eukaryotes all stay represented)
    recs = [ctx.rec(k) for k in LOW_CARRIER if ctx.rec(k).usable_for_positives]
    nid = _ids("M-planted")
    i = 0
    for t_sub in (1, 2, 3):
        for m in range(t_sub, 4):
            for j in range(3):
                # the shortest length whose ball keeps the expected random contexts low enough for the list to be
                # released whole (the smoke run: a 20-mer at max_mismatches 2 has 21,563 level-2 contexts), with the
                # window's composition (_planted_window, G03)
                rec = recs[(i * 7 + j * 3 + 1) % len(recs)]
                i += 1
                tr = _planted_window(ctx, rng, rec, [[20, 25, 31][(i - 1 + j) % 3], 25, 31], m,
                                     ctx.rng("M", "plant-redraw", i))
                L = tr["length"]
                r2 = ctx.rng("M", "plant", i)
                positions = st.spread_positions(r2, L, t_sub)
                variant, edits = st.substitute_bases(tr["window"], positions, r2)
                tr["edits"] = edits
                tr["variant"] = variant
                bits = st.mismatch_information_bits_iupac(variant, m)
                if bits < 24:
                    continue
                req = pattern_request([_p(variant, "v")], max_mismatches=m)
                subs = _subs_plus(variant, tr["window"], positions)
                out.append(scenario(nid(), "M", "%d-mer of %s with %d planted substitution(s), max_mismatches %d: the source at level %d" % (L, rec.key, t_sub, m, t_sub),
                                    "POST /pattern", req, tr,
                                    [OK, INV, CC, ex("contains_source_context", strand="+", require="all", level=t_sub, substitutions=subs),
                                     ex("count_relation", count="contexts", level=t_sub, relation="exact", min=len(source_contexts(tr)))],
                                    tags=["mismatch", "planted", "machine-dependent"]
                                    + (["smoke"] if (t_sub == 2 and m == 2 and j == 0) else [])))
    # the reverse orientation: P = rc(variant), the source a - context with substitutions in P's coordinates
    for j in range(4):
        rec = recs[(j * 5 + 4) % len(recs)]
        tr = _planted_window(ctx, rng, rec, [25, 31], 2, ctx.rng("M", "rc-redraw", j))
        Lr = tr["length"]
        r2 = ctx.rng("M", "rc", j)
        positions = st.spread_positions(r2, Lr, 2)
        variant, edits = st.substitute_bases(tr["window"], positions, r2)
        P = st.revcomp(variant)
        tr["edits"] = edits
        tr["strand"] = "-"
        req = pattern_request([_p(P, "rc")], max_mismatches=2)
        out.append(scenario(nid(), "M", "reverse orientation: the source a - context at level 2, substitutions mirrored", "POST /pattern",
                            req, tr, [OK, INV, CC, ex("contains_source_context", strand="-", require="all", level=2,
                                                                  substitutions=_subs_minus(P, tr["window"]))],
                            tags=["mismatch", "planted", "reverse", "machine-dependent"]))
    # m < t: the source is out of reach (its level is t)
    nid = _ids("M-below")
    for j in range(6):
        rec = recs[(j * 9 + 2) % len(recs)]
        tr = _mm_window(ctx, rng, rec, 31)
        r2 = ctx.rng("M", "below", j)
        t_sub = 2 + j % 2
        positions = st.spread_positions(r2, 31, t_sub)
        variant, edits = st.substitute_bases(tr["window"], positions, r2)
        tr["edits"] = edits
        req = pattern_request([_p(variant, "v")], mode="partial", max_mismatches=t_sub - 1)
        out.append(scenario(nid(), "M", "%d planted substitutions, max_mismatches %d: the source not at any level" % (t_sub, t_sub - 1),
                            "POST /pattern", req, tr, [OK, INV, CC, ex("source_absent", strand="+")],
                            tags=["mismatch", "near-miss"]))
    # IUPAC + substitutions
    nid = _ids("M-iupac")
    for j in range(4):
        rec = recs[(j * 7 + 6) % len(recs)]
        tr = _mm_window(ctx, rng, rec, 25)
        r2 = ctx.rng("M", "iupac", j)
        w = list(tr["window"])
        ipos = sorted(r2.sample(range(25), 3))
        for p in ipos:
            w[p] = st.iupac_code_holding(w[p], 2, r2) if p != ipos[-1] else "N"
        free = [p for p in range(25) if p not in ipos]
        spos = sorted(r2.sample(free, 1))
        for p in spos:
            w[p] = r2.choice([b for b in "ACGT" if b != tr["window"][p]])
        pat = "".join(w)
        if st.mismatch_information_bits_iupac(pat, 1) < 24:
            continue
        tr["iupac_positions"] = ipos
        tr["edits"] = [{"pos": p, "from": tr["window"][p], "to": pat[p]} for p in spos]
        subs = [{"position": p, "expected": pat[p], "found": tr["window"][p]} for p in spos]
        req = pattern_request([_p(pat, "iu", kind="iupac")], max_mismatches=1)
        out.append(scenario(nid(), "M", "IUPAC codes + one substitution: level 1, never at N, a base inside a set no substitution",
                            "POST /pattern", req, tr,
                            [OK, INV, CC, ex("contains_source_context", strand="+", require="all", level=1, substitutions=subs)],
                            tags=["mismatch", "iupac"]))
    # protein: one substituted residue (a sense codon of another residue), m = 1
    from truth import load_proteins
    # machine-dependent: the smoke validation of lane S1b saw a table-4 decapeptide at max_mismatches 1 stopped by its
    # 60 s budget in discovery on the loaded host (withheld deadline): the scenario is then inconclusive (v6 gives
    # every default-budget /pattern request 110,000 ms, PATTERN_TIME_BUDGET_MS; M-protein-0010 keeps its own 30,000)
    nid = _ids("M-protein")
    prots = [p for p in load_proteins() if p["genome"] in ctx.by_key and p["transl_table"] in (1, 11)]
    n = 0
    for prot in prots[::7]:
        win = next((w for w in prot["windows"] if w["length_aa"] == 10), None)
        if win is None:
            continue
        tr = _peptide_truth(ctx, prot, win)
        if tr is None:
            continue
        r2 = ctx.rng("M", "prot", prot["key"])
        res_i = r2.randrange(1, 9)
        orig = win["peptide"][res_i]
        new = r2.choice([a for a in st.AMINO if a != orig])
        pv = win["peptide"][:res_i] + new + win["peptide"][res_i + 1:]
        codon = win["dna"][3 * res_i:3 * res_i + 3]
        tr["residue_edit"] = {"residue": res_i, "from": orig, "to": new, "record_codon": codon}
        subs = [{"residue": res_i, "expected": new, "found": orig, "codon": codon}]
        req = pattern_request([{"id": "pv", "protein": pv}], genetic_code=prot["transl_table"], max_mismatches=1)
        out.append(scenario(nid(), "M", "decapeptide of %s with residue %d %s->%s, max_mismatches 1: the source at level 1" % (prot["key"], res_i, orig, new),
                            "POST /pattern", req, tr,
                            [OK, INV, CC, ex("contains_source_context", strand=tr["strand"], require="all", level=1, substitutions=subs)],
                            timing_class="slow", tags=["mismatch", "protein", "machine-dependent"]))
        n += 1
        if n >= 6:
            break
    # protein positives in a code other than 11 (the verification's M4), and at m = 2 under today's cap of 3
    np_ = 0
    nc_tables: List[int] = []
    for prot in [p for p in load_proteins() if p["genome"] in ctx.by_key and p["transl_table"] in (4, 12)]:
        win = next((w for w in prot["windows"] if w["length_aa"] == 10 and st.translate(w["dna"], 1) != w["peptide"]), None)
        tr = _peptide_truth(ctx, prot, win) if win else None
        if tr is None:
            continue
        r2 = ctx.rng("M", "prot-nc", prot["key"])
        res_i = r2.randrange(1, 9)
        orig = win["peptide"][res_i]
        new = r2.choice([a for a in st.AMINO if a != orig])
        pv = win["peptide"][:res_i] + new + win["peptide"][res_i + 1:]
        codon = win["dna"][3 * res_i:3 * res_i + 3]
        tr["residue_edit"] = {"residue": res_i, "from": orig, "to": new, "record_codon": codon}
        subs = [{"residue": res_i, "expected": new, "found": orig, "codon": codon}]
        req = pattern_request([{"id": "pv", "protein": pv}], genetic_code=prot["transl_table"], max_mismatches=1)
        out.append(scenario(nid(), "M", "decapeptide of %s (table %d) with residue %d %s->%s, max_mismatches 1: the source "
                            "at level 1" % (prot["key"], prot["transl_table"], res_i, orig, new), "POST /pattern", req, tr,
                            [OK, INV, CC, ex("contains_source_context", strand=tr["strand"], require="all",
                                                         level=1, substitutions=subs)],
                            timing_class="slow", tags=["mismatch", "protein", "machine-dependent",
                                                       "table-%d" % prot["transl_table"]]))
        np_ += 1
        nc_tables.append(prot["transl_table"])
        if np_ >= 2:
            break
    if 12 not in nc_tables:
        ctx.skip("M", "M-protein table 12", "no table-12 decapeptide whose coding DNA reads differently in table 1 "
                 "(table 12's proteins hold no CTG in any window or shift); the scenarios built: table %s" % (
                     ", ".join(str(x) for x in nc_tables) or "none"))
    # m = 2 only under a protein cap of at least 2 in the build's capabilities snapshot (G04: the deployed cap is 1,
    # and above it the host answers the cap slot, the next scenario); the ids are kept either way
    pcap = HOST_MISMATCH_MAX.get("protein", 3)
    for j, prot in enumerate(prots[1::11][:2]):
        if pcap < 2:
            ctx.skip("M", "%s (m = 2 decapeptide of %s)" % (nid(), prot["key"]),
                     "max_mismatches 2 is above the snapshot's protein cap of %d (its cap slot is the next scenario)" % pcap)
            continue
        win = next((w for w in prot["windows"] if w["length_aa"] == 10), None)
        tr = _peptide_truth(ctx, prot, win) if win else None
        if tr is None:
            ctx.skip("M", "M-protein m=2 %s" % prot["key"], "its decapeptide is not in a fetched record")
            continue
        r2 = ctx.rng("M", "prot2", prot["key"])
        res = sorted(r2.sample(range(1, 9), 2))
        pv = list(win["peptide"])
        subs = []
        for res_i in res:
            orig = pv[res_i]
            pv[res_i] = r2.choice([a for a in st.AMINO if a != orig])
            subs.append({"residue": res_i, "expected": pv[res_i], "found": orig, "codon": win["dna"][3 * res_i:3 * res_i + 3]})
        tr["residue_edits"] = subs
        req = pattern_request([{"id": "pv2", "protein": "".join(pv)}], genetic_code=prot["transl_table"], max_mismatches=2)
        out.append(scenario(nid(), "M", "decapeptide of %s with 2 substituted residues, max_mismatches 2 (the snapshot's "
                            "protein cap %d): the source at level 2" % (prot["key"], pcap), "POST /pattern", req, tr,
                            [OK, INV, CC, ex("contains_source_context", strand=tr["strand"], require="all",
                                                         level=2, substitutions=subs)],
                            timing_class="slow", tags=["mismatch", "protein", "machine-dependent"]))
    # suffix scope with substitutions: one offset, the floor without exemption (§21.5)
    nid_s = _ids("M-suffix")
    for j, key in enumerate(("hvo", "pfa")):
        rec = ctx.rec(key)
        tr = _mm_window(ctx, rng, rec, 18)
        r2 = ctx.rng("M", "suffix", j)
        positions = st.spread_positions(r2, 18, 1)
        variant, edits = st.substitute_bases(tr["window"], positions, r2)
        tr["edits"], tr["variant"] = edits, variant
        if st.mismatch_information_bits_iupac(variant, 1) < 24:
            ctx.skip("M", "M-planted suffix %s" % key, "below the floor at m 1")
            continue
        req = pattern_request([_p(variant, "v")], scope="suffix", max_mismatches=1)
        out.append(scenario(nid_s(), "M", "an 18-mer of %s with 1 substitution in suffix scope, max_mismatches 1: the source "
                            "at offset k-L, level 1" % key, "POST /pattern", req, tr,
                            [OK, INV, CC, ex("contains_source_context", strand="+", require="all", level=1,
                                             substitutions=_subs_plus(variant, tr["window"], positions))],
                            tags=["mismatch", "planted", "suffix"]))
    # the snapshot's protein cap below 3: max_mismatches cap + 1 for a peptide -> its slot (the service refuses the
    # submit: mismatch.above_kind_cap, G05); above a cap of 3 a value is refused whatever the kind (M-refusal)
    if prots:
        sid = nid()
        if pcap >= 3:
            ctx.skip("M", "%s (a peptide above the protein cap)" % sid,
                     "the snapshot's protein cap is %d: no max_mismatches above it is in range" % pcap)
        else:
            prot = prots[0]
            win = next(w for w in prot["windows"] if w["length_aa"] == 10)
            tr = _peptide_truth(ctx, prot, win)
            mm = pcap + 1
            req = pattern_request([{"id": "cap", "protein": win["peptide"]}], genetic_code=prot["transl_table"],
                                  max_mismatches=mm, max_contexts=PATTERN_MAX_CONTEXTS, time_budget_ms=30000)
            out.append(scenario(sid, "M", "a peptide with max_mismatches %d above the protein cap of %d" % (mm, pcap),
                                "POST /pattern", req, tr,
                                [OK, ex("pattern_slot_error", pattern=0, code="mismatches_above_cap",
                                        text="pattern: max_mismatches %d above this server's cap of %d for protein patterns "
                                             "(the capabilities' mismatches.max.protein); ask with at most %d" % (mm, pcap, pcap),
                                        fields={"max_mismatches_cap": pcap})],
                                timing_class="slow", tags=["mismatch", "protein"]))
    # order and stops
    nid = _ids("M-order")
    rec = ctx.rec("eco")
    tr = _mm_window(ctx, rng, rec, 20)
    req = pattern_request([_p(tr["window"])], mode="partial", max_contexts=5, max_mismatches=1)
    # the record alone holds 12 level-0 contexts of its 20-mer, so a list cut at 5 holds level 0 only (mismatches
    # first); which 5 of them, in node order, the truth cannot say (the verification's S2)
    out.append(scenario(nid(), "M", "a 20-mer, max_mismatches 1, partial max_contexts 5: the 5 released all at level 0",
                        "POST /pattern", req, tr,
                        [OK, INV, CC, ex("results_levels", max_level=0, min_returned=5), ex("cut", reason="max_contexts"),
                         ex("count_relation", level=0, relation="exact", min=len(source_contexts(tr)))],
                        tags=["mismatch", "order"]))
    tr = _mm_window(ctx, rng, rec, 31)
    req = pattern_request([_p(tr["window"])], mode="count", max_mismatches=2, max_steps=1500)
    out.append(scenario(nid(), "M", "max_steps 1500 at max_mismatches 2: stop in a level, lower levels exact", "POST /pattern",
                        req, tr, [OK, INV, CC, ex("stop", phase="discovery", reason="max_steps", mismatch_level_in=[1, 2]),
                                  ex("count_relation", level=0, relation="exact", min=1)], tags=["mismatch", "stop"]))
    req = pattern_request([_p(tr["window"])], mode="partial", max_mismatches=2, max_steps=1500)
    out.append(scenario(nid(), "M", "the same stop in partial: the lower levels released, cut max_steps", "POST /pattern",
                        req, tr, [OK, INV, CC, ex("cut", reason="max_steps"),
                                  ex("contains_source_context", strand="+", require="all", level=0, substitutions=[])],
                        tags=["mismatch", "stop"]))
    # the floor at m
    nid = _ids("M-floor")
    for L, m in ((17, 2), (14, 1), (19, 3)):
        tr = _mm_window(ctx, rng, rec, L)
        bits = st.mismatch_information_bits_iupac(tr["window"], m)
        within = None
        for mm in range(m - 1, -1, -1):
            ok = (st.mismatch_information_bits_iupac(tr["window"], mm) >= 24) if mm else (2 * L >= 24)
            if ok:
                within = mm
                break
        text = "pattern: %.1f information bits within %d substitution%s, below the floor of 24.0 for scope any_offset" % (
            bits, m, "" if m == 1 else "s")
        req = pattern_request([_p(tr["window"])], mode="count", max_mismatches=m)
        out.append(scenario(nid(), "M", "exact %d-mer at max_mismatches %d: below the floor, %s within it" % (L, m, within),
                            "POST /pattern", req, tr,
                            [OK, ex("pattern_slot_error", pattern=0, code="information_below_floor", text=text,
                                    fields={"max_mismatches_within_floor": within})],
                            tags=["mismatch", "floor"] + (["smoke"] if L == 17 else [])))
    # refusals: the host's 400 text; the MCP coerces or refuses these values (async: not comparable through it)
    nid = _ids("M-refusal")
    tr = _mm_window(ctx, rng, rec, 25)
    txt = "request.max_mismatches: expected an integer from 0 to 3 (substitutions per instance)"
    for v in (4, -1, 1.5, "1", None, True):
        req = pattern_request([_p(tr["window"])], max_mismatches=None)
        req["max_mismatches"] = v
        out.append(scenario(nid(), "M", "max_mismatches %s: refused" % repr(v), "POST /pattern", req, tr,
                            [err400(txt, code="invalid_request")], tags=["mismatch", "refusal"]))
    later = ("request.max_mismatches: substitutions with output.labels \"all\" or \"predicate_only\", a predicate or "
             "long_search \"supported_paths\" are not served by this build (the capabilities' mismatches.later_increment); "
             "ask without them, or with max_mismatches 0")
    for extra, title in (({"output": {"labels": "all"}}, "labels all"),
                         ({"predicate": {"any": [rec.column]}}, "a predicate"),
                         ({"long_search": "supported_paths"}, "supported paths")):
        req = pattern_request([_p(tr["window"])], max_mismatches=1, **extra)
        out.append(scenario(nid(), "M", "max_mismatches 1 with %s: later_increment" % title, "POST /pattern", req, tr,
                            [err400(later, code="later_increment")], tags=["mismatch", "refusal"]))
    t40 = _mm_window(ctx, rng, rec, 40)
    for ls in (None, "paths"):
        req = pattern_request([_p(t40["window"], "long"), _p(tr["window"], "short")], mode="count", max_mismatches=2,
                              long_search=ls)
        out.append(scenario(nid(), "M", "a 40-mer at max_mismatches 2 (long_search %s): its slot, cap 0" % ls, "POST /pattern",
                            req, t40,
                            [OK, ex("pattern_slot_error", pattern=0, code="mismatches_above_cap",
                                    text="pattern: max_mismatches 2 for a pattern of 40 bases, longer than k (31): substitutions in patterns longer than k are not served by this build (the capabilities' mismatches.later_increment names long_patterns); ask with max_mismatches 0, or with windows of at most k bases",
                                    fields={"max_mismatches_cap": 0}),
                             ex("count_relation", pattern=1, level=0, relation="exact", min=7)],
                            tags=["mismatch", "slots"]))
    req = pattern_request([_p(tr["window"])], max_mismatches=0)
    out.append(scenario(nid(), "M", "max_mismatches 0: the answer without the field", "POST /pattern", req, tr,
                        [OK, INV, CC, ex("contains_source_context", strand="+", require="all"),
                         ex("answer_fields", absent=["limits.max_mismatches"]),
                         ex("entry_fields", absent=["counts.contexts.by_mismatches", "max_mismatches"])],
                        tags=["mismatch"]))
    return out


# ================================================================================================= Q: predicates
def gen_Q(ctx: Ctx) -> List[dict]:
    out: List[dict] = []
    rng = ctx.rng("Q")
    nid = _ids("Q-pred")
    human = ctx.rec("hsa").column
    for j, key in enumerate(LOW_CARRIER[:6]):
        rec = ctx.rec(key)
        t = ctx.window(rng, rec, 25)
        A = rec.column
        cases = [
            ({"any": [A]}, True, "any(A) selects the record's contexts"),
            ({"none": [A]}, False, "none(A) deselects them"),
            ({"all": [A, human]}, False, "all(A, human) on a window human lacks"),
            ({"at_least": {"n": 1, "labels": [A, human]}}, True, "at_least 1 of (A, human)"),
            ({"at_least": {"n": 2, "labels": [A, human]}}, False, "at_least 2 of (A, human)"),
            ({"and": [{"any": [A]}, {"not": {"any": [human]}}]}, True, "and(any A, not any human)"),
            ({"or": [{"any": [human]}, {"any": [A]}]}, True, "or(any human, any A)"),
        ]
        pred, sel, title = cases[j % len(cases)]
        req = pattern_request([_p(t["window"])], output={"labels": "predicate_only"}, predicate=pred,
                              predicate_strands="context")
        exp = [OK, INV, ex("predicate_selected", strand="+", selected=sel)]
        if sel:
            exp.append(ex("predicate_selected", strand="+", selected=True, selection_strand_in=["context", "both"], label=A))
        out.append(scenario(nid(), "Q", "%s on %s" % (title, key), "POST /pattern", req, t, exp,
                            tags=["predicate"] + (["smoke"] if j == 1 else [])))
    for j, key in enumerate(LOW_CARRIER[6:10]):
        rec = ctx.rec(key)
        t = ctx.window(rng, rec, 25)
        A = rec.column
        pred, sel, title = [({"any": [A]}, True, "any(A)"), ({"none": [A]}, False, "none(A)"),
                            ({"not": {"any": [A]}}, False, "not(any A)"), ({"all": [A]}, True, "all(A)")][j]
        req = pattern_request([_p(t["window"])], output={"labels": "predicate_only"}, predicate=pred,
                              predicate_strands="either")
        out.append(scenario(nid(), "Q", "%s under predicate_strands either on %s" % (title, key), "POST /pattern", req, t,
                            [OK, INV, ex("predicate_selected", strand="+", selected=sel)], tags=["predicate"]))
    rec = ctx.rec("buc")
    t = ctx.window(rng, rec, 25)
    req = pattern_request([_p(t["window"])], predicate={"any": ["999999999"]})
    out.append(scenario(nid(), "Q", "an unknown label: listed, any() folds to false, nothing read", "POST /pattern", req, t,
                        [OK, INV, ex("answer_fields", equals={"predicate.unknown_labels": ["999999999"], "predicate.normal_form": False}),
                         ex("entry_fields", equals={"selection.pass": "constant", "counts.selected.value": 0}),
                         ex("notes", has=["predicate_constant"])], tags=["predicate", "unknown-label"]))
    req = pattern_request([_p(t["window"])], predicate={"none": ["999999999"]})
    out.append(scenario(nid(), "Q", "none(unknown) folds to true: every context selected", "POST /pattern", req, t,
                        [OK, INV, ex("answer_fields", equals={"predicate.normal_form": True}),
                         ex("predicate_selected", strand="+", selected=True)], tags=["predicate", "unknown-label"]))
    req = pattern_request([_p(t["window"])], predicate={"any": [562]})
    out.append(scenario(nid(), "Q", "a number as a name: refused with the fix", "POST /pattern", req, t,
                        [err400("request.predicate.any[0]: expected a string (a column label; write a taxid as \"562\")",
                                code="invalid_request")], tags=["predicate", "refusal"]))
    req = pattern_request([_p(t["window"])], output={"labels": "predicate_only"})
    out.append(scenario(nid(), "Q", "predicate_only without a predicate: refused", "POST /pattern", req, t,
                        [err400("request.output.labels: \"predicate_only\" returns the labels a predicate names: it needs a predicate (request.predicate)",
                                code="invalid_request")], tags=["predicate", "refusal"]))
    # Q-strands: the mirror under "either", on 20-mers of a 16S gene whose 42 bp around them the record also holds
    # in reverse (rRNA operons on both strands): the '-' contexts are then the reverse complements of the record's
    # own k-mers, which must exist (the verification's A2: a mirror check with no mirror tests nothing)
    nid = _ids("Q-strands")
    for j, key in enumerate(["eco", "sau", "pae"]):
        rec = ctx.rec(key)
        t = _inverted_window(ctx, rec, 20)
        if t is None:
            ctx.skip("Q", "Q-strands %s" % key, "no 16S window with an inverted copy of its 42 bp")
            continue
        req = pattern_request([_p(t["window"])], output={"labels": "predicate_only"}, predicate={"any": [rec.column]},
                              predicate_strands="either")
        out.append(scenario(nid(), "Q", "either on a 16S 20-mer of %s held on both strands: the + contexts by their own "
                            "row, the mirrors (which exist) by their reverse complement" % key,
                            "POST /pattern", req, t,
                            [OK, INV, ex("predicate_selected", strand="+", selected=True, selection_strand_in=["context", "both"], label=rec.column),
                             ex("predicate_either_mirror", label=rec.column, min=1)],
                            timing_class="slow", tags=["predicate", "strands"]))
    # a window taken from the reverse strand under predicate_strands context: its source contexts are '-' and
    # selected by their own rows; `or` and an at_least above the list's size under either
    nid = _ids("Q-more")
    rq = ctx.rng("Q", "more")
    rec = ctx.rec("mja")
    t = ctx.window(rq, rec, 25)
    t["strand"] = "-"
    req = pattern_request([_p(st.revcomp(t["window"]))], output={"labels": "predicate_only"},
                          predicate={"any": [rec.column]}, predicate_strands="context")
    out.append(scenario(nid(), "Q", "the reverse complement of a mja window, context strands: the '-' source contexts "
                        "selected by their own rows", "POST /pattern", req, t,
                        [OK, INV, ex("predicate_selected", strand="-", selected=True, selection_strand_in=["context"],
                                     label=rec.column)], tags=["predicate", "strands"]))
    rec = ctx.rec("tth")
    t = ctx.window(rq, rec, 25)
    req = pattern_request([_p(t["window"])], output={"labels": "predicate_only"},
                          predicate={"or": [{"any": [human]}, {"any": [rec.column]}]}, predicate_strands="either")
    out.append(scenario(nid(), "Q", "or(any human, any A) under either on tth: selected", "POST /pattern", req, t,
                        [OK, INV, ex("predicate_selected", strand="+", selected=True)], tags=["predicate"]))
    rec = ctx.rec("hvo")
    t = ctx.window(rq, rec, 25)
    # n above the list's length is a form fault: SPEC-pattern-search §19.3 (1 <= n <= the list's length), refused
    # naming the path of the first fault on both arms (G07; the service: pattern.predicate_form_other)
    req = pattern_request([_p(t["window"])], output={"labels": "predicate_only"},
                          predicate={"at_least": {"n": 3, "labels": [rec.column, human]}}, predicate_strands="either")
    out.append(scenario(nid(), "Q", "at_least 3 of two labels: refused, n outside 1 .. the list's length",
                        "POST /pattern", req, t,
                        [err400(code="invalid_request", contains=["request.predicate.at_least.n"])],
                        tags=["predicate", "refusal"]))
    # predicates on supported paths (Q-supported): none() deselects the record's path, any() under either selects it
    nid = _ids("Q-supported")
    for j, key in enumerate(("buc", "sso")):
        rec = ctx.rec(key)
        t = ctx.window(rq, rec, 60)
        req = pattern_request([_p(t["window"])], mode="partial", long_search="supported_paths",
                              output={"labels": "predicate_only"}, predicate={"none": [rec.column]},
                              predicate_strands="context", max_paths=1024)
        out.append(scenario(nid(), "Q", "supported paths of a %s window under none(A): the record's path not selected" % key,
                            "POST /pattern", req, t, [OK, INV, ex("path_absent", sequence=t["window"])],
                            tags=["predicate", "supported-paths"]))
        req = pattern_request([_p(t["window"])], long_search="supported_paths", output={"labels": "predicate_only"},
                              predicate={"any": [rec.column]}, predicate_strands="either", max_paths=1024)
        out.append(scenario(nid(), "Q", "supported paths of the same window under any(A), either strands: the record "
                            "verified", "POST /pattern", req, t, [OK, INV, ex("path_record_verified", strand="+")],
                            tags=["predicate", "supported-paths", "strands"]))
    # Q-motif
    nid = _ids("Q-motif")
    rec = ctx.rec("mge")
    t = ctx.window(rng, rec, 25)
    req = pattern_request([_p(t["window"])], mode="count", predicate={"any": [rec.column]}, predicate_scope="motif")
    out.append(scenario(nid(), "Q", "motif any(A) on a window of A: selected, every_context", "POST /pattern", req, t,
                        [OK, INV, ex("motif", selected=True, decided_by_in=["every_context"], label_present=rec.column)],
                        tags=["predicate", "motif", "smoke"]))
    q = st.random_dna(ctx.rng("Q", "motif-absent"), 31)
    req = pattern_request([_p(q)], mode="count", predicate={"any": [rec.column]}, predicate_scope="motif")
    out.append(scenario(nid(), "Q", "motif on an absent pattern: no_instance", "POST /pattern", req, {"random": q},
                        [OK, INV, ex("motif", selected=False, decided_by_in=["no_instance"])], tags=["predicate", "motif", "negative"]))
    req = pattern_request([_p(t["window"])], mode="count", predicate={"none": [rec.column]}, predicate_scope="motif")
    out.append(scenario(nid(), "Q", "motif none(A) on a window of A: false, every_context", "POST /pattern", req, t,
                        [OK, INV, ex("motif", selected=False, decided_by_in=["every_context"])], tags=["predicate", "motif"]))
    req = pattern_request([_p(t["window"])], mode="count", predicate_scope="motif")
    out.append(scenario(nid(), "Q", "predicate_scope motif without a predicate: refused", "POST /pattern", req, t,
                        [err400(code="invalid_request", contains=["predicate_scope"])], tags=["predicate", "refusal"]))
    # predicate_strands without a predicate: a no-op on the host, a refusal of the service (A-12 b)
    req = pattern_request([_p(t["window"])], predicate_strands="context")
    out.append(scenario(nid(), "Q", "predicate_strands without a predicate: a no-op, the source found", "POST /pattern",
                        req, t, [OK, INV, CC, ex("contains_source_context", strand="+", require="all")],
                        tags=["predicate", "split-case"]))
    # Q-budget
    nid = _ids("Q-budget")
    rec = ctx.rec("cje")
    t = ctx.window(rng, rec, 25)
    req = pattern_request([_p(t["window"])], predicate={"any": [rec.column]}, max_predicate_contexts=1)
    out.append(scenario(nid(), "Q", "max_predicate_contexts 1: not admitted", "POST /pattern", req, t,
                        [OK, INV, ex("withheld", reason="predicate_above_threshold")], tags=["predicate", "budget"]))
    req = pattern_request([_p(t["window"])], mode="partial", predicate={"any": [rec.column]}, max_predicate_contexts=1)
    out.append(scenario(nid(), "Q", "max_predicate_contexts 1 in partial: the first tested, cut", "POST /pattern", req, t,
                        [OK, INV, ex("cut", reason="max_predicate_contexts")], tags=["predicate", "budget"]))
    req = pattern_request([_p(t["window"])], predicate={"any": [rec.column]}, max_predicate_work=1, predicate_strands="context")
    out.append(scenario(nid(), "Q", "max_predicate_work 1: the selection stops", "POST /pattern", req, t,
                        [OK, INV, ex("withheld", reason="predicate_budget"),
                         ex("stop", phase="selection", reason="max_predicate_work")], tags=["predicate", "budget"]))
    return out


_16S_MOTIF = "GTGCCAGCAGCCGCGGTAA"   # the 515F site, exact in the 16S genes of eco, sau and pae


def _inverted_window(ctx: Ctx, rec, L: int) -> Optional[dict]:
    """A window of L bases downstream of a forward 16S site whose surrounding (k - L) bases on each side, every
    k-mer holding the window, also occur reverse-complemented in the record (an operon on the other strand)."""
    seq = rec.seq
    pad = K - L
    for p in st.find_all(seq, _16S_MOTIF):
        for off in range(30, 900, 37):
            s0 = p + off
            region = seq[s0 - pad:s0 + L + pad]
            if len(region) == L + 2 * pad and st.is_acgt(region) and st.revcomp(region) in seq \
                    and st.low_complexity_dust_score(seq[s0:s0 + L]) < 2:
                t = ctx.truth_at(rec, s0, L)
                t["inverted_copy"] = rec.offset + seq.find(st.revcomp(region))
                return t
    return None


# ================================================================================================= C: capabilities
def gen_C(ctx: Ctx, caps: Dict[str, dict]) -> List[dict]:
    out: List[dict] = []
    nid = _ids("C-caps")
    tc = caps["traverse_capabilities"]
    pc = caps["pattern_capabilities"]
    out.append(scenario(nid(), "C", "GET /traverse/capabilities: identity, limits, the gate block within 32 KiB",
                        "GET /traverse/capabilities", None, None,
                        [OK, ex("capabilities_doc", required=["release", "index_fp", "k", "max_query_bp", "max_seeds",
                                                              "pattern.pattern_contract_version", "pattern.details",
                                                              "attempts.server_instance"],
                                equals={"index_fp": tc["index_fp"], "k": 31, "regime": "basic",
                                        "pattern.pattern_contract_version": 1, "pattern.available": True,
                                        "pattern.details": "GET /pattern/capabilities"},
                                max_bytes=32768)], tags=["capabilities", "smoke"],
                        # the MCP's capabilities tool answers the service's view (gate, backends, service_maxima,
                        # presets), not the host's document: not implemented via MCP (the owner, 2026-10-15)
                        mcp={"not_implemented": "the host's /traverse/capabilities document (the MCP capabilities "
                                                "tool answers the service's view, not the host's document)",
                             "by_design": ["capabilities.traverse_capabilities_document"]}))
    out.append(scenario(nid(), "C", "GET /pattern/capabilities: the full block", "GET /pattern/capabilities", None, None,
                        [OK, ex("capabilities_doc", required=["caps", "mismatches", "predicate", "genetic_codes", "default_max_anchors"],
                                equals={"pattern_contract_version": 1, "available": True, "graph_mode": "basic", "k": 31,
                                        "mask": "file", "counting": "exact", "placement": "record",
                                        "support": "record_verified", "default_projection": "none"})],
                        tags=["capabilities", "smoke"]))
    out.append(scenario(nid(), "C", "GET /capabilities: features, routes, the pattern block, ready", "GET /capabilities", None, None,
                        [OK, ex("capabilities_doc", required=["features", "routes", "pattern", "release", "server_instance", "resolve"],
                                equals={"ready": True, "mode": "single", "routes.pattern": "POST /pattern",
                                        "routes.pattern_capabilities": "GET /pattern/capabilities"})],
                        tags=["capabilities"]))
    nid = _ids("C-protocol")
    eco = ctx.rec("eco")
    t = ctx.window(ctx.rng("C"), eco, 25)
    body = pattern_request([_p(t["window"])])
    out.append(scenario(nid(), "C", "Accept-Encoding gzip: a gzip answer", "POST /pattern", body, t,
                        [OK, ex("content_encoding", value="gzip"), ex("contains_source_context", strand="+", require="all")],
                        tags=["protocol"], headers={"Accept-Encoding": "gzip"}))
    out.append(scenario(nid(), "C", "Accept-Encoding deflate: a deflate answer", "POST /pattern", body, t,
                        [OK, ex("content_encoding", value="deflate")], tags=["protocol"], headers={"Accept-Encoding": "deflate"}))
    out.append(scenario(nid(), "C", "an error body is never compressed", "POST /pattern", dict(body, bogus=1), t,
                        [err400("request: unknown field 'bogus'", code="invalid_request"),
                         ex("content_encoding", value=None, any_status=True)],
                        tags=["protocol"], headers={"Accept-Encoding": "gzip"}))
    out.append(scenario(nid(), "C", "a JSON text with a comment: refused", "POST /pattern", None, None,
                        [err400(code="invalid_request")], tags=["protocol", "refusal"],
                        body_raw='{"patterns": [{"dna": "%s"}] /* x */}' % t["window"]))
    out.append(scenario(nid(), "C", "a duplicated member name: refused", "POST /pattern", None, None,
                        [err400(code="invalid_request")], tags=["protocol", "refusal"],
                        body_raw='{"patterns": [{"dna": "%s"}], "mode": "count", "mode": "count"}' % t["window"]))
    out.append(scenario(nid(), "C", "a body that is not an object: refused", "POST /pattern", [1, 2], None,
                        [err400(code="invalid_request")], tags=["protocol", "refusal"]))
    out.append(scenario(nid(), "C", "a wrong content type: the JSON body is still read", "POST /pattern", body, t,
                        [OK, ex("contains_source_context", strand="+", require="all")], tags=["protocol", "code-derived"],
                        headers={"Content-Type": "text/plain"}))
    out.append(scenario(nid(), "C", "GET on POST /pattern: not served", "GET /pattern", None, None,
                        [ex("status", http=[400, 404, 405])], tags=["protocol"]))
    out.append(scenario(nid(), "C", "GET on POST /resolve: not served", "GET /resolve", None, None,
                        [ex("status", http=[400, 404, 405])], tags=["protocol"]))
    return out


# ================================================================================================= D: determinism
def _status_200_only(s: dict) -> bool:
    """Whether the scenario's status expectation admits a 200 only (manifest_lib.allowed_statuses)."""
    return allowed_statuses(s.get("expect")) == [200]


def gen_D(ctx: Ctx, pool: List[dict]) -> List[dict]:
    """Twenty fast scenarios of other areas, repeated three times: identical answers apart from timing."""
    import copy as _copy
    out: List[dict] = []
    picks = []
    want = {"P": 6, "R": 4, "L": 3, "M": 3, "T": 2, "Q": 2}
    for area, n in want.items():
        # a twin repeats its source's request three times: never one whose answer depends on what ran before it (an
        # attempt id, run-unique and retained: a D twin would meet its source's 409, P4-builder-truth-03), one that
        # repeats already, or one whose status allows another answer than a 200
        cands = [s for s in pool if s["area"] == area and s["timing_class"] == "fast" and "refusal" not in s["tags"]
                 and "machine-dependent" not in s["tags"] and "retention-dependent" not in s["tags"]
                 and "relation" not in s["tags"] and "split-case" not in s["tags"] and "warm-up" not in s["tags"]
                 and "gate-coverage" not in s["tags"] and not s.get("run_unique") and s.get("repeat", 1) == 1
                 and _status_200_only(s) and not s.get("body_raw") and not s.get("headers")]
        step = max(1, len(cands) // n)
        picks += cands[::step][:n]
    # a complete search at m >= 2 (A-12 a: the machine-dependent tag keeps planted searches out of the picks above):
    # the forward planted window at max_mismatches 2 with the lowest chance estimate, so it completes fast (a deadline
    # repeat reads inconclusive, G27)
    m2 = [s for s in pool if s["id"].startswith("M-planted") and "reverse" not in s["tags"]
          and (s.get("request") or {}).get("max_mismatches") == 2]
    if m2:
        picks.append(min(m2, key=lambda s: ((s.get("truth") or {}).get("chance_contexts_estimate", 1e18), s["id"])))
    for j, s in enumerate(picks):
        d = _copy.deepcopy(s)
        d["id"] = "D-%s" % s["id"]
        d["area"] = "D"
        d["title"] = "determinism: %s" % s["title"]
        d["repeat"] = 3
        d["expect"] = [OK, ex("deterministic")] + [e for e in s["expect"] if e["check"] != "status"]
        d["tags"] = sorted(set(s["tags"]) - {"smoke"} | {"determinism"} | ({"smoke"} if j in (0, 6) else set()))
        d["derived_from"] = s["id"]
        out.append(d)
    return out


# ================================================================================================= N: negatives, near misses, hard regions, use cases
PRIMERS = {   # universal 16S primers (IUPAC), the classic agent use case
    "27F": "AGAGTTTGATCMTGGCTCAG", "1492R": "GGTTACCTTGTTACGACTT",
    "515F": "GTGYCAGCMGCCGCGGTAA", "806R": "GGACTACNVGGGTWTCTAAT",
}


def _iupac_sites(seq: str, pat: str) -> List[Tuple[int, str]]:
    """Every (start, strand) of the IUPAC pattern on both strands of |seq| (brute force over a few Mb)."""
    L = len(pat)
    out = []
    rc = st.revcomp(pat)
    anchor = pat[:6] if all(c in "ACGT" for c in pat[:6]) else None
    for s in range(0, len(seq) - L + 1):
        w = seq[s:s + L]
        if st.iupac_match(pat, w):
            out.append((s, "+"))
        if st.iupac_match(rc, w):
            out.append((s, "-"))
    return out


def _tandem_repeats(seq: str, min_len: int = 40) -> List[Tuple[int, int, str]]:
    """Perfect tandem repeats of period 1-6 (homopolymers included) of at least |min_len| bases in |seq| (ACGT
    units): (start, length, unit), the longest first, each once at its smallest period."""
    found = []
    for p in range(1, 7):
        run, j, n = 0, p, len(seq)
        while j <= n:
            if j < n and seq[j] == seq[j - p]:
                run += 1
            else:
                if run + p >= min_len:
                    start = j - run - p
                    unit = seq[start:start + p]
                    if st.is_acgt(unit) and len(set(unit)) > (1 if p > 1 else 0):
                        found.append((start, run + p, unit))
                run = 0
            j += 1
    # keep a repeat once, at its smallest period
    found.sort(key=lambda x: (-x[1], len(x[2])))
    out, seen = [], set()
    for st0, ln, unit in found:
        if st0 not in seen:
            out.append((st0, ln, unit))
            seen.add(st0)
    return out


def _real_repeats(ctx: Ctx) -> List[dict]:
    """Low-complexity windows found in the fetched records (the verification's M8: DESIGN N asks for real ones):
    the longest perfect tandem repeat of each record that has one of at least 40 bp, a 31-mer inside it as a
    pattern (the record's contexts, a count at least theirs, the note where sdust flags it), 60 bp of it resolved
    against the record's column, and two traversed under the record's header (a walk in a repeat re-enters its
    own edges: walk_in_record's exact minimum states where)."""
    out: List[dict] = []
    nid = _ids("N-realrep")
    walked = 0
    for rec in ctx.positives:
        reps = [r for r in _tandem_repeats(rec.seq) if r[0] > 100 and r[0] + 100 < len(rec.seq)]
        if not reps:
            continue
        s0, ln, unit = reps[0]
        t = ctx.truth_at(rec, s0 + 2, K)
        t["tandem_repeat"] = {"start": rec.offset + s0, "length": ln, "unit": unit}
        lc = st.sdust_flags(t["window"])
        out.append(scenario(nid(), "N", "a 31-mer inside %s's %d bp (%s)n repeat: the record's contexts, the note %s" % (
                                rec.key, ln, unit, "stated" if lc else "absent"), "POST /pattern",
                            pattern_request([_p(t["window"], "rep")]), t,
                            [OK, INV, CC, ex("contains_source_context", strand="+", require="all"),
                             ex("count_relation", relation="exact", min=len(set(source_contexts(t)))),
                             ex("low_complexity_note", expected=lc)],
                            tags=["hard-region", "low-complexity", "real-repeat"]))
        if ln >= 60:
            t60 = ctx.truth_at(rec, s0, 60)
            out.append(scenario(nid(), "N", "60 bp of %s's (%s)n repeat resolved: its column supports every k-mer" % (
                                    rec.key, unit), "POST /resolve", resolve_request(t60["window"], labels=[rec.column]),
                                t60, [OK, INV, ex("resolve_runs", label=rec.column, full=True)],
                                tags=["hard-region", "real-repeat"]))
        if walked < 2 and not rec.region and ln >= 60:
            ts = ctx.flank_for(rec, ctx.truth_at(rec, s0, 40), 250)
            exp = [OK, INV, ex("walk_in_record", label_scope="record")]
            if len(unit) <= 40 - K:
                # a period of at most 40 - k bases repeats a seed k-mer one period on inside the 40 bp seed, so each
                # arm's step into the repeat's next copy is a step into a seed node (rejoined_seed at 0 bp, SPEC
                # §6.6): the re-entry the title claims, which a walk minimum of 0 bp leaves untested
                # (P5-v4-content-02)
                exp.append(ex("blocked_event", reason="rejoined_seed", max_at_bp=0))
            out.append(scenario(nid(), "N", "a 40 bp seed in %s's (%s)n repeat walked under its header: the record "
                                "followed to where it re-enters or branches" % (rec.key, unit), "POST /traverse",
                                traverse_request([_tseed(ts["window"], [rec.accession])], {"bounds": {"max_extension_bp": 200}}),
                                ts, exp, tags=["hard-region", "real-repeat", "traverse"]))
            walked += 1
    return out


def gen_N(ctx: Ctx) -> List[dict]:
    return _within_kind_caps(ctx, "N", _gen_N(ctx))


def _gen_N(ctx: Ctx) -> List[dict]:
    out: List[dict] = []
    rng = ctx.rng("N")
    recs = ctx.positives
    eco = ctx.rec("eco")
    # negatives on every route
    nid = _ids("N-negative")
    for j, L in enumerate((31, 40, 60, 100, 200)):
        q = st.random_dna(ctx.rng("N", "rand", L), L)
        if L <= K:
            req = pattern_request([_p(q)])
            exp = [OK, INV, CC, ex("absent")]
        else:
            req = pattern_request([_p(q)], long_search="paths")
            exp = [OK, INV, CC, ex("absent", count="paths")]
        out.append(scenario(nid(), "N", "random %d bp on /pattern: exact 0, complete" % L, "POST /pattern", req, {"random": q},
                            exp, tags=["negative"]))
    for j in range(4):
        rec = recs[(j * 7 + 3) % len(recs)]
        t = ctx.window(rng, rec, 60)
        q = st.shuffled(ctx.rng("N", "shuf", j), t["window"])
        out.append(scenario(nid(), "N", "shuffled 60 bp of %s on /pattern: exact 0 paths" % rec.key, "POST /pattern",
                            pattern_request([_p(q)], long_search="paths"), {"shuffled_from": t},
                            [OK, INV, CC, ex("absent", count="paths")], tags=["negative", "shuffled"] + (["smoke"] if j == 0 else [])))
        out.append(scenario(nid(), "N", "shuffled 60 bp of %s on /resolve: no k-mer in the graph" % rec.key, "POST /resolve",
                            resolve_request(q, labels=[rec.column]), {"shuffled_from": t},
                            [OK, INV, ex("resolve_runs", label=rec.column, kmers=0)], tags=["negative", "shuffled"]))
    for j, L in enumerate((100, 200)):
        q = st.random_dna(ctx.rng("N", "rres", L), L)
        out.append(scenario(nid(), "N", "random %d bp on /resolve, discovery: nothing to discover" % L, "POST /resolve",
                            resolve_request(q, discover={"max_labels": 20, "kind": "column"}), {"random": q},
                            [OK, INV, ex("resolve_absent")], tags=["negative", "discover"]))
        out.append(scenario(nid(), "N", "random %d bp seed on /traverse: the seed's 400, never a walk" % L, "POST /traverse",
                            traverse_request([_tseed(q)], {"bounds": {"max_extension_bp": 100}}), {"random": q},
                            [err400(prefix="seed '%s': Seed is not fully present in the graph; graph runs: " % q[:32])],
                            tags=["negative", "refusal"]))
    # the shortest a random peptide of these residues is whose chance bound stays below 1e-5 is 12 residues (36 bases,
    # paths; 10 residues, a pattern below k, met about 5.6e-4: chance_bound)
    for j, n_aa in enumerate((12, 20, 30)):
        r2 = ctx.rng("N", "pep", n_aa)
        pep = "".join(r2.choice("MWCFYHQNKDE") for _ in range(n_aa))
        if 3 * n_aa <= K:
            req = pattern_request([{"id": "rp", "protein": pep}], genetic_code=1)
            exp = [OK, INV, CC, ex("absent")]
        else:
            req = pattern_request([{"id": "rp", "protein": pep}], genetic_code=1, long_search="paths")
            exp = [OK, INV, CC, ex("absent", count="paths")]
        out.append(scenario(nid(), "N", "a random %d-residue peptide of low-degeneracy residues: exact 0" % n_aa, "POST /pattern",
                            req, {"random_peptide": pep}, exp, tags=["negative", "peptide"]))
    # near misses: t = m+1 .. m+5 substitutions, spread or clustered; pairs reached / out of reach
    nid = _ids("N-nearmiss")
    for j in range(16):
        rec = recs[(j * 5 + 1) % len(recs)]
        t = ctx.window(rng, rec, 31)
        r2 = ctx.rng("N", "near", j)
        m = j % 4
        t_sub = m + 1 + (j // 4) % 5
        t_sub = max(4, min(t_sub, 8)) if j % 2 == 0 else m + 1 + (j % 5)
        t_sub = min(t_sub, 8)
        clustered = j % 3 == 0
        positions = st.clustered_positions(r2, 31, t_sub) if clustered else st.spread_positions(r2, 31, t_sub)
        variant, edits = st.substitute_bases(t["window"], positions, r2)
        t["edits"] = edits
        t["variant"] = variant
        mm = {} if m == 0 else {"max_mismatches": m}
        req = pattern_request([_p(variant, "near")], mode="partial", **mm)
        out.append(scenario(nid(), "N", "%d %s substitutions in a 31-mer of %s, max_mismatches %d: the source out of reach" % (
            t_sub, "clustered" if clustered else "spread", rec.key, m), "POST /pattern", req, t,
            [OK, INV, CC, ex("source_absent", strand="+")],
            tags=["near-miss", "machine-dependent"] + (["smoke"] if j == 5 else [])))
        if m and t_sub > m:
            # its pair within reach: the same window with m of the substitutions
            keep = sorted(positions[:m])
            v2, e2 = st.substitute_bases(t["window"], [], r2)
            v2 = list(t["window"])
            for e in edits:
                if e["pos"] in keep:
                    v2[e["pos"]] = e["to"]
            v2 = "".join(v2)
            t2 = dict(t)
            t2["edits"] = [e for e in edits if e["pos"] in keep]
            t2["variant"] = v2
            subs = _subs_plus(v2, t["window"], keep)
            req = pattern_request([_p(v2, "reach")], mode="partial", max_mismatches=m)
            out.append(scenario(nid(), "N", "its pair with %d substitution(s): the source back at level %d" % (m, m), "POST /pattern",
                                req, t2, [OK, INV, CC, ex("contains_source_context", strand="+", require="all",
                                                                      level=m, substitutions=subs)],
                                tags=["near-miss", "pair", "machine-dependent"]))
    # peptides with 2-4 substituted residues under the protein cap (1 after the next restart; 3 today)
    from truth import load_proteins
    prots = [p for p in load_proteins() if p["genome"] in ctx.by_key and p["transl_table"] in (1, 11)]
    for j, prot in enumerate(prots[3::17][:4]):
        win = next((w for w in prot["windows"] if w["length_aa"] == 10), None)
        tr = _peptide_truth(ctx, prot, win) if win else None
        if tr is None:
            continue
        r2 = ctx.rng("N", "pepnear", j)
        k_sub = 2 + j % 3
        res = sorted(r2.sample(range(10), k_sub))
        pv = list(win["peptide"])
        for r in res:
            pv[r] = r2.choice([a for a in st.AMINO if a != pv[r]])
        pv = "".join(pv)
        tr["residue_edits"] = res
        req = pattern_request([{"id": "pn", "protein": pv}], genetic_code=prot["transl_table"], mode="partial", max_mismatches=1)
        out.append(scenario(nid(), "N", "decapeptide with %d substituted residues at max_mismatches 1: the source out of reach" % k_sub,
                            "POST /pattern", req, tr, [OK, INV, CC, ex("source_absent", strand=tr["strand"])],
                            timing_class="slow", tags=["near-miss", "peptide", "machine-dependent"]))
    # partly mutated seeds for /resolve and /traverse: a substitution every ~k/2 bases
    for j in range(3):
        rec = recs[(j * 11 + 4) % len(recs)]
        t = ctx.window(rng, rec, 200)
        w = list(t["window"])
        sub_pos = list(range(15, 200, 15))
        for p in sub_pos:
            w[p] = "A" if w[p] != "A" else "C"
        q = "".join(w)
        t["edits"] = [{"pos": p} for p in sub_pos]
        out.append(scenario(nid(), "N", "200 bp of %s with a substitution every 15 bases: no k-mer of the record survives" % rec.key,
                            "POST /resolve", resolve_request(q, labels=[rec.column]), t,
                            [OK, INV, ex("resolve_runs", label=rec.column, max_kmers=0)], tags=["near-miss", "mutated-seed"]))
        q2 = t["window"][:100] + ("A" if t["window"][100] != "A" else "C") + t["window"][101:]
        # the 31 k-mers holding the substitution: on a whole record its header supports exactly the k-mers the record
        # holds (exact runs), while the taxid column may hold the variant in another record (pae's column holds 1,762:
        # P5-v4-content-03), so the column is asked to cover both sides only
        labels = [rec.column] + ([] if rec.region else [rec.accession])
        exp = [OK, INV, ex("resolve_runs", label=rec.column, covers=[0, 70]),
               ex("resolve_runs", label=rec.column, covers=[101, 170])]
        if not rec.region:
            exp.append(ex("resolve_runs", label=rec.accession, runs=_record_runs(rec, q2)))
        out.append(scenario(nid(), "N", "200 bp with one substitution at 100: the runs break there, both sides supported",
                            "POST /resolve", resolve_request(q2, labels=labels), t, exp, tags=["near-miss", "mutated-seed"]))
        # the seed under the record's header: a seed refusal, either one (as T-strand): "not fully present" when a
        # variant k-mer is in no record of the graph, "no seed label supports every k-mer" when another record holds
        # them all; the header supports none of them where the whole record holds none (checked here)
        held = [i for i in range(70, 101) if q2[i:i + K] in rec.seq]
        tt = dict(t, variant_kmers_in_record=len(held),
                  seed_refusals=["seed '%s': Seed is not fully present in the graph; graph runs: " % q2[:32],
                                 "seed '%s': No seed label supports every k-mer of the seed" % q2[:32]])
        tags = ["near-miss", "mutated-seed", "refusal"]
        if rec.region:
            tt = empirical(tt, "the record holds none of the variant k-mers outside its fetched region")
            tags.append("index-empirical")
        if held:
            ctx.skip("N", "a traverse seed with a substitution in %s" % rec.key,
                     "the record holds %d of the variant's k-mers elsewhere: its header would support the seed" % len(held))
        else:
            out.append(scenario(nid(), "N", "a traverse seed with a substitution in its middle, under the record's "
                                "header: a seed refusal",
                                "POST /traverse", traverse_request([_tseed(q2, [rec.accession])], {"bounds": {"max_extension_bp": 100}}),
                                tt, [err400(prefix="seed '%s': " % q2[:32])], tags=tags))
    # hard regions: homopolymers, dinucleotide and short tandem repeats, microsatellites, multi-copy, telomeres.
    # Contexts are graph k-mers, not occurrences: a pattern of L <= k has at most 2 (k - L + 1) 4^(k - L) of them
    # (A^31: 2), so the counts are exact and bounded, the patterns occur in refseq33m, and the low_complexity_pattern
    # note is asserted exactly where sdust flags the pattern (seqtools.sdust_flags, the host's T = 20, W = 64)
    nid = _ids("N-hard")
    hard = [("A" * 31, "A^31"), ("A" * 40, "A^40"), ("AT" * 15 + "A", "(AT)^15A"), ("AT" * 16, "(AT)^16"),
            ("CAG" * 10, "(CAG)^10"), ("CA" * 12 + "C", "(CA)^12C"), ("TTAGGG" * 5, "(TTAGGG)^5 telomeric"),
            ("CCCTAA" * 5, "(CCCTAA)^5 telomeric")]
    for q, title in hard:
        lc = st.sdust_flags(q)
        note = ex("low_complexity_note", expected=lc)
        truth = {"low_complexity": title, "sdust_flags": lc}
        if len(q) > K:
            req = pattern_request([_p(q, "hard")], mode="partial", long_search="paths", max_paths=100)
            out.append(scenario(nid(), "N", "%s (%d bp > k) with paths, partial: its path counted, the note %s" % (
                                    title, len(q), "stated" if lc else "absent"), "POST /pattern", req, truth,
                                [OK, INV, CC, ex("count_relation", count="paths", relation="exact", min=1), note],
                                tags=["hard-region", "low-complexity", "paths"]))
            continue
        bound = 2 * (K - len(q) + 1) * 4 ** (K - len(q))
        req = pattern_request([_p(q, "hard")], mode="partial", max_contexts=100)
        out.append(scenario(nid(), "N", "%s (%d bp) in partial: exact, at most %d contexts, the note %s" % (
                                title, len(q), bound, "stated" if lc else "absent"), "POST /pattern", req, truth,
                            [OK, INV, CC, ex("count_relation", relation="exact", min=1, max=bound), note],
                            tags=["hard-region", "low-complexity"] + (["smoke"] if title == "(AT)^15A" else [])))
        req = pattern_request([_p(q, "hard")], max_contexts=100)
        if bound <= 100:
            # at most `bound` contexts, never above max_contexts 100: the list is complete, never withheld
            out.append(scenario(nid(), "N", "%s in all_or_count with max_contexts 100: at most %d contexts, so complete "
                                "and never withheld, exact" % (title, bound), "POST /pattern", req, truth,
                                [OK, INV, CC, ex("withheld", reason=None),
                                 ex("count_relation", relation="exact", min=1, max=bound), note],
                                tags=["hard-region", "low-complexity"]))
            continue
        out.append(scenario(nid(), "N", "%s in all_or_count with max_contexts 100: complete (or withheld above 100), "
                            "exact" % title, "POST /pattern", req, truth,
                            [OK, INV, CC, ex("withheld", reason=[None, "count_above_threshold"]),
                             ex("count_relation", relation="exact", min=1, max=bound), note],
                            tags=["hard-region", "low-complexity"]))
    for q, title in (("A" * 80, "A^80"), ("TTAGGG" * 15, "(TTAGGG)^15"), ("CAG" * 30, "(CAG)^30")):
        req = pattern_request([_p(q, "hardlong")], mode="partial", long_search="paths", max_paths=100, max_anchors=1000)
        lc = st.sdust_flags(q)
        exp = [OK, INV, CC, ex("low_complexity_note", expected=lc)]
        if q == "A" * 80:
            exp.append(ex("count_relation", count="paths", relation="exact", min=1))
        out.append(scenario(nid(), "N", "%s with paths: stated stops or cut, the note %s" % (title, "stated" if lc else "absent"),
                            "POST /pattern", req, {"low_complexity": title, "sdust_flags": lc}, exp,
                            tags=["hard-region", "low-complexity", "paths"]))
    for j, (key, t, copies) in enumerate(_repeat_windows(ctx, ctx.rng("N", "rep"), n=2)):
        rec = ctx.rec(key)
        w31 = t["window"][:31]
        tt = ctx.truth_at(rec, t["start"] - rec.offset, 31)
        req = pattern_request([_p(w31, "rep")], mode="count")
        out.append(scenario(nid(), "N", "a 31-mer of a %d-copy repeat of %s: counted" % (len(copies), key), "POST /pattern", req, tt,
                            [OK, INV, ex("count_relation", relation="exact", min=1)], tags=["hard-region", "repeat"]))
        out.append(scenario(nid(), "N", "resolve a %d-copy repeat of %s: every k-mer carried" % (len(copies), key), "POST /resolve",
                            resolve_request(t["window"], labels=[rec.column]), t,
                            [OK, INV, ex("resolve_runs", label=rec.column, full=True)], tags=["hard-region", "repeat"]))
    q = "TTAGGG" * 10
    # the seed's period closes a cycle inside the seed (its 30 k-mers use the cycle's 6 edges), so each arm's first
    # step into the next copy is a step into a seed node: a blocked successor `rejoined_seed` at 0 bp on both arms
    # (SPEC-labeled-traversal-core §6.6), whatever else the human column holds around telomeric repeats. Its 6
    # distinct k-mers are in nearly every vertebrate record: their rows took the cold MCP call past its attempt's
    # bound in the seed phase (17.5 s, triage-v4 G24), so a warm-up resolve of the seed under the human column (the
    # rows the seed phase reads) runs right before it
    hard_id = nid()
    hsa_col = ctx.rec("hsa").column
    out.append(scenario(hard_id, "N", "a telomeric seed (TTAGGG)^10 traversed: each arm re-enters the seed at 0 bp "
                        "(rejoined_seed), with its stated outcome (or the seed's 400 if the human column lacks it)",
                        "POST /traverse",
                        traverse_request([_tseed(q, [hsa_col])], {"bounds": {"max_extension_bp": 200, "time_budget_ms": 5000}}),
                        {"low_complexity": "telomeric", "period": 6, "copies_in_seed": 10},
                        [ex("status", http=[200, 400]), INV, ex("traverse_seed_outcome", walks=["complete", "partial"]),
                         ex("blocked_event", reason="rejoined_seed", max_at_bp=0)],
                        tags=["hard-region", "telomeric", "traverse"]))
    hard_warm = warm_up(ctx, WARM_IDS[hard_id], q, hsa_col,
                        empirical({"low_complexity": "telomeric", "period": 6},
                                  "the human column carries every k-mer of (TTAGGG)^10 (v3 and v4: the seed walked on "
                                  "every arm)"),
                        hard_id, "the rows of the telomeric repeat's k-mers in the human column",
                        tags=("index-empirical", "telomeric"))
    out += _real_repeats(ctx)
    # every normal use case of an agent, with its failure modes
    nid = _ids("N-usecase")
    sites = {}
    for name, primer in PRIMERS.items():
        sites[name] = _iupac_sites(eco.seq, primer)
    for name, primer in PRIMERS.items():
        ss = sites[name]
        if not ss:
            continue
        s, strand = ss[0]
        tt = ctx.truth_at(eco, s, len(primer))
        tt["strand"] = strand
        n_src = len({(km, p) for (s2, strand2) in ss for (km, p) in source_contexts(ctx.truth_at(eco, s2, len(primer)))})
        req = pattern_request([_p(primer, name, kind="iupac")], mode="count")
        out.append(scenario(nid(), "N", "primer %s (IUPAC) on both strands: counted, at least E. coli K-12's %d site contexts" % (name, n_src),
                            "POST /pattern", req, tt,
                            [OK, INV, ex("count_relation", relation="exact", min=n_src)],
                            tags=["use-case", "primer"] + (["smoke"] if name == "515F" else [])))
    req = pattern_request([_p(PRIMERS["515F"], "515F", kind="iupac"), _p(PRIMERS["806R"], "806R", kind="iupac")],
                          mode="count", strands="forward")
    # both entries searched forward only and both counted (P5-v4-content-04: the second primer was not asserted, so its
    # slot refused would have passed); 806R's forward count is the records that hold the 16S gene on the reverse
    # strand (v4: 6,900 exact, 515F 7,737)
    out.append(scenario(nid(), "N", "a primer pair on the forward strand only (one primer is reverse): both counted", "POST /pattern",
                        req, None, [OK, INV] + [e for j in (0, 1) for e in (
                            ex("entry_fields", pattern=j, equals={"strands": ["+"]}),
                            ex("count_relation", pattern=j, relation="exact", min=1))],
                        tags=["use-case", "primer", "wrong-strand"]))
    t = ctx.window(rng, eco, 25)
    req = pattern_request([_p(st.revcomp(t["window"]), "rc")], strands="forward")
    out.append(scenario(nid(), "N", "the wrong strand: the reverse complement of a window searched forward only: the record's k-mers not found",
                        "POST /pattern", req, dict(t, strand="+"),
                        [OK, INV, CC, ex("source_absent", strand="+")], tags=["use-case", "wrong-strand"]))
    req = pattern_request([_p(t["window"][:10], "short")], mode="count")
    out.append(scenario(nid(), "N", "too short (10 bp, 20 bits): below the floor", "POST /pattern", req, t,
                        [OK, ex("pattern_slot_error", pattern=0, code="information_below_floor",
                                text="pattern: 20.0 information bits, below the floor of 24.0 for scope any_offset")],
                        tags=["use-case", "too-short"]))
    req = pattern_request([_p("N" * 12 + t["window"][12:25], "degenerate", kind="iupac")], mode="count")
    out.append(scenario(nid(), "N", "degenerate (12 N + 13 bases, 26 bits): admitted; the 13 bases at every offset are "
                        "costly (200 with its counts, or 503 deadline on a loaded host)", "POST /pattern", req, t,
                        [OK, INV, ex("count_relation", relation=["exact", "at_least"], min=1)],
                        timing_class="slow", tags=["use-case", "degenerate", "machine-dependent", "heavy"]))
    req = pattern_request([_p("NNNNNNNNNNNNNNNNNNNN" + t["window"][:11], "degenerate2", kind="iupac")], mode="count")
    out.append(scenario(nid(), "N", "too degenerate (20 N + 11 bases, 22 bits): below the floor", "POST /pattern", req, t,
                        [OK, ex("pattern_slot_error", pattern=0, code="information_below_floor",
                                text="pattern: 22.0 information bits, below the floor of 24.0 for scope any_offset")],
                        tags=["use-case", "degenerate"]))
    # identify an unknown read: resolve + pattern of the same window
    for j in range(2):
        rec = recs[(j * 13 + 8) % len(recs)]
        t = ctx.window(rng, rec, 150)
        out.append(scenario(nid(), "N", "identify a read of %s: discovery names its taxid" % rec.key, "POST /resolve",
                            resolve_request(t["window"], discover={"max_labels": 100, "kind": "column"}), t,
                            [OK, INV, ex("resolve_discover_contains", label=rec.column)], timing_class="slow", tags=["use-case", "identify"]))
        t31 = ctx.truth_at(rec, t["start"] - rec.offset, 31)
        out.append(scenario(nid(), "N", "identify a read of %s: its first 31 bp placed in its record" % rec.key, "POST /pattern",
                            pattern_request([_p(t31["window"])], output={"labels": "predicate_only"},
                                            predicate={"any": [rec.column]}, predicate_strands="context"), t31,
                            [OK, INV, ex("source_placement", strand="+")], tags=["use-case", "identify"]))
    # a peptide in the wrong genetic code (mitochondrial table 2 vs table 1)
    from truth import load_proteins as _lp
    mito = [p for p in _lp() if p["transl_table"] in (2, 3, 4, 5) and p["genome"] in ctx.by_key]
    for prot in mito[:2]:
        win = next((w for w in prot["windows"] if w["length_aa"] == 10 and st.translate(w["dna"], 1) != w["peptide"]), None)
        tr = _peptide_truth(ctx, prot, win) if win else None
        if tr is None:
            # stated, not dropped (P5-v4-content-05); the X-table1 scenarios test the other codes' wrong-code reads
            ctx.skip("N", "N-usecase own table / table 1 pair of %s (table %d)" % (prot["key"], prot["transl_table"]),
                     "no decapeptide of this protein reads differently in table 1" if win is None else
                     "its differing decapeptide is not in a fetched record")
            continue
        out.append(scenario(nid(), "N", "a %s peptide in its own table %d: found" % (prot["genome"], prot["transl_table"]), "POST /pattern",
                            pattern_request([{"id": "own", "protein": win["peptide"]}], genetic_code=prot["transl_table"]), tr,
                            [OK, INV, CC, ex("contains_source_context", strand=tr["strand"], require="all")], tags=["use-case", "genetic-code"]))
        out.append(scenario(nid(), "N", "the same peptide in table 1: its record codons out of reach", "POST /pattern",
                            pattern_request([{"id": "wrong", "protein": win["peptide"]}], genetic_code=1), tr,
                            [OK, INV, CC, ex("source_absent", strand=tr["strand"])], tags=["use-case", "genetic-code", "wrong-code"]))
    # extend a seed to its neighbourhood (both arms) and a long marker as supported paths: in T and L
    # the index defect, stated: the A. thaliana chloroplast (last record of column 3702) is not annotated
    athc = ctx.rec("athc")
    t = ctx.truth_at(athc, 54989, 31)
    out.append(scenario(nid(), "N", "the 2026-10-09 defect: rbcL 31-mer of NC_000932.1 is in the graph but its header carries none",
                        "POST /resolve", resolve_request(t["window"], labels=[athc.accession]), t,
                        [OK, INV, ex("resolve_graph_runs", runs=[[0, 1]]), ex("resolve_runs", label=athc.accession, kmers=0)],
                        tags=["index-defect", "documented"]))
    out += _near_miss_more(ctx)
    _annotate_walks(ctx, out)
    out.append(hard_warm)   # placed right before its walk by the build (build_manifest.place_warm_ups)
    return out


def _near_miss_more(ctx: Ctx) -> List[dict]:
    """Near misses at another length and on the minus strand (the verification's M8), on low-carrier records: a
    25-mer with 3 substitutions at max_mismatches 1 (out of reach), the reverse complement of a 31-mer with 4 at 2
    (out of reach, the source a '-' context), and its pair with 1 at 1 (the source back at level 1, the
    substitutions in P's coordinates, SPEC §21.8)."""
    out: List[dict] = []
    rng = ctx.rng("N", "near-more")
    nid = _ids("N-nearmiss2")
    rec = ctx.rec("cal")
    t = ctx.window(rng, rec, 25)
    pos = st.spread_positions(rng, 25, 3)
    v, edits = st.substitute_bases(t["window"], pos, rng)
    t["edits"], t["variant"] = edits, v
    out.append(scenario(nid(), "N", "a 25-mer of cal with 3 spread substitutions, max_mismatches 1: the source out of reach",
                        "POST /pattern", pattern_request([_p(v, "near")], mode="partial", max_mismatches=1), t,
                        [OK, INV, CC, ex("source_absent", strand="+")], tags=["near-miss", "machine-dependent"]))
    rec = ctx.rec("mja")
    t = ctx.window(rng, rec, 31)
    pos = st.spread_positions(rng, 31, 4)
    v, edits = st.substitute_bases(t["window"], pos, rng)
    tm = dict(t, strand="-", edits=edits, variant=v)
    out.append(scenario(nid(), "N", "the reverse complement of a mja 31-mer with 4 substitutions, max_mismatches 2: the "
                        "'-' source out of reach", "POST /pattern",
                        pattern_request([_p(st.revcomp(v), "near")], mode="partial", max_mismatches=2), tm,
                        [OK, INV, CC, ex("source_absent", strand="-")],
                        tags=["near-miss", "reverse", "machine-dependent"]))
    keep = [edits[0]]
    v1 = list(t["window"])
    v1[keep[0]["pos"]] = keep[0]["to"]
    v1 = "".join(v1)
    P = st.revcomp(v1)
    t1 = dict(t, strand="-", edits=keep, variant=v1)
    out.append(scenario(nid(), "N", "its pair with 1 substitution at max_mismatches 1: the '-' source back at level 1",
                        "POST /pattern", pattern_request([_p(P, "reach")], mode="partial", max_mismatches=1), t1,
                        [OK, INV, CC, ex("contains_source_context", strand="-", require="all", level=1,
                                                     substitutions=_subs_minus(P, t["window"]))],
                        tags=["near-miss", "pair", "reverse", "machine-dependent"]))
    return out


# ================================================================================================= chance bounds
# A global negative control (a random or shuffled sequence, a random peptide) claims an absence the records cannot
# certify: the chance that the graph holds the query anyway is bounded, and the builder refuses a control whose bound
# exceeds CHANCE_BOUND_MAX (VERIFY-D of the outside review)
CHANCE_BOUND_MAX = 1e-4
_NEGATIVE_TRUTH_KEYS = ("random", "shuffled_from", "random_peptide")


def chance_bound(sc: dict) -> Optional[dict]:
    """The chance bound of a global negative control (None for another scenario), the union bound of VERIFY-D over the
    graph's REFSEQ33M_KMERS k-mers (k = 31) for a uniform query: a query k-mer equals some graph k-mer with probability
    at most N / 4^k. /pattern, DNA of L >= k: no context or path spells the query on a requested strand, which needs its
    first k-mer: (strands) N / 4^k; a peptide of L = 3n bases over its E encodings in its table: below k a substring of
    some k-mer, (strands) E (k - L + 1) N / 4^L, from k on a path of L bases, of which the graph spells about as many
    as it holds k-mers (a model, not a strict bound): (strands) E N / 4^L. /resolve: any of the query's k-mers in the
    graph, (L - k + 1) N / 4^k (a column's support bounded by the whole graph). /traverse (the seed's 400): every
    k-mer present, bounded by the first, N / 4^k. A shuffled query also states the composition-matched worst case
    (the whole graph at the query's base composition: the bound times (4 sum f^2)^k)."""
    t = sc.get("truth")
    req = sc.get("request")
    if not isinstance(t, dict) or not any(k in t for k in _NEGATIVE_TRUTH_KEYS) or not isinstance(req, dict):
        return None
    n = REFSEQ33M_KMERS
    pk = n / 4.0 ** K
    route = sc.get("route")
    out: Dict[str, Any] = {"kmers": n, "k": K}
    q = None
    if route == "POST /pattern":
        p0 = (req.get("patterns") or [{}])[0]
        ns = 1 if req.get("strands") in ("forward", "reverse") else 2
        if "protein" in p0:
            pep, table = p0["protein"], int(req.get("genetic_code", 1))
            enc = 1
            for a in pep:
                enc *= len(st.codons_of(a, table))
            L = 3 * len(pep)
            p = ns * enc * ((K - L + 1) * n / 4.0 ** L if L < K else n / 4.0 ** L)
            out.update(model="peptide of %d bases, %d encodings in table %d, %d strand(s)" % (L, enc, table, ns))
        else:
            q = p0.get("dna") or p0.get("iupac") or ""
            L = len(q)
            p = ns * ((K - L + 1) * n / 4.0 ** L if L < K else pk)
            out.update(model="%d bases, %d strand(s): the first k-mer" % (L, ns))
    elif route == "POST /resolve":
        q = req.get("sequence") or ""
        p = (len(q) - K + 1) * pk
        out.update(model="any of %d k-mers in the graph" % (len(q) - K + 1))
    elif route == "POST /traverse":
        q = ((req.get("seeds") or [{}])[0]).get("sequence") or ""
        p = pk
        out.update(model="every k-mer of the seed in the graph, bounded by the first")
    else:
        return None
    out["p"] = float("%.3g" % p)
    if "shuffled_from" in t and q:
        f = [q.count(b) / float(len(q)) for b in "ACGT"]
        worst = (4.0 * sum(x * x for x in f)) ** K
        out["p_composition_worst"] = float("%.3g" % min(1.0, p * worst))
    return out


def state_chance_bounds(scenarios: List[dict]) -> None:
    """truth.chance_bound on every global negative control (chance_bound); a truth shared by two scenarios is copied
    first, so that each states its own route's bound."""
    for sc in scenarios:
        b = chance_bound(sc)
        if b is not None:
            sc["truth"] = dict(sc["truth"], chance_bound=b)


# ================================================================================================= relations
def _by_title(scenarios: List[dict], area: str, title: str) -> dict:
    hits = [x for x in scenarios if x["area"] == area and x["title"] == title]
    if len(hits) != 1:
        raise ValueError("%d %s scenarios titled %r" % (len(hits), area, title))
    return hits[0]


def _add_tag(sc: dict, tag: str) -> None:
    sc["tags"] = sorted(set(sc["tags"]) | {tag})


def _by_id(scenarios: List[dict], sid: str) -> dict:
    return next(x for x in scenarios if x["id"] == sid)


def _twin(src: dict, sid: str, title: str, request, expect: List[dict], truth: Optional[dict] = None,
          tags=(), timing_class: Optional[str] = None) -> dict:
    """A scenario that differs from |src| only in what the relation compares (tagged relation)."""
    return scenario(sid, src["area"], title, src["route"], request, src["truth"] if truth is None else truth, expect,
                    timing_class=timing_class or src["timing_class"],
                    tags=sorted(set(tags) | {"relation"}))


def gen_REL(ctx: Ctx, scenarios: List[dict], caps: Dict[str, dict]) -> List[dict]:
    """The manifest's relations (DESIGN §3: R-budget's prefix re-resolve and its continuation, suffix within
    any_offset, an IUPAC pattern's contexts holding the exact window's, max_mismatches 0 answered as without the
    field, the three capabilities documents consistent) and the twin scenarios they need. Runs after every other
    area but D (whose picks leave relation twins out), on random streams of its own, so that the windows of the
    other scenarios stay as built. Relations are recorded on ctx.relations; the twins are returned."""
    out: List[dict] = []
    # same_answer: max_mismatches 0 is answered as the request without the field (SPEC §21.2)
    m0 = _by_title(scenarios, "M", "max_mismatches 0: the answer without the field")
    req = copy.deepcopy(m0["request"])
    del req["max_mismatches"]
    tw = _twin(m0, "M-twin-0001", "the m=0 window without max_mismatches", req,
               [OK, INV, CC, ex("contains_source_context", strand="+", require="all")], tags=["mismatch", "smoke"])
    _add_tag(m0, "smoke")
    out.append(tw)
    relation(ctx, "same_answer", m0["id"], tw["id"], "max_mismatches 0 answers byte for byte as the request without it")
    # subset: the exact answer's contexts are the m=1 answer's level-0 contexts
    req = copy.deepcopy(m0["request"])
    req["max_mismatches"] = 1
    tw1 = _twin(m0, "M-twin-0002", "the m=0 window at max_mismatches 1: its level 0 holds the exact answer", req,
                [OK, INV, CC, ex("contains_source_context", strand="+", require="all", level=0, substitutions=[])],
                tags=["mismatch"])
    out.append(tw1)
    relation(ctx, "subset", tw["id"], tw1["id"], "the exact contexts are among the level-0 contexts at max_mismatches 1",
             a_pattern=0, b_pattern=0, key="context", b_level=0)
    # same_answer: lower case is upper-cased (pattern and resolve)
    lc = _by_title(scenarios, "P", "lower case: upper-cased, the record found")
    req = copy.deepcopy(lc["request"])
    for p in req["patterns"]:
        p["dna"] = p["dna"].upper()
    tw = _twin(lc, "P-twin-0001", "the lower-case window in upper case", req,
               [OK, INV, ex("contains_source_context", strand="+", require="all")])
    out.append(tw)
    relation(ctx, "same_answer", lc["id"], tw["id"], "a lower-case pattern answers as its upper case")
    lc = _by_title(scenarios, "R", "lower case: resolved as upper case")
    req = copy.deepcopy(lc["request"])
    req["sequence"] = req["sequence"].upper()
    lab = req["labels"][0]
    tw = _twin(lc, "R-twin-0001", "the lower-case resolve in upper case", req,
               [OK, INV, ex("resolve_runs", label=lab, full=True)])
    out.append(tw)
    relation(ctx, "same_answer", lc["id"], tw["id"], "a lower-case resolve answers as its upper case")
    # prefix_equals: a budget stop answers the resolve of a prefix; the continuation from a k-mer completes it. The
    # complete twin runs right before its budget scenario (a warm-up: build_manifest.place_warm_ups) under the default
    # 60,000 ms, which the cold query met within (v4: 15.6 s for eco's second half, 21.9 s for pae's 10 kb)
    for i, sid in enumerate(("R-budget-0001", "R-budget-0002")):
        rb = _by_id(scenarios, sid)
        seq = rb["request"]["sequence"]
        lab = rb["request"]["labels"][0]
        base = {k: v for k, v in rb["truth"].items() if k != "stop_design"}
        full = _twin(rb, "R-twin-%04d" % (2 + 2 * i), "the same %s bp under the %s ms budget, run first (it warms the "
                     "rows): the profile the stop must be a prefix of" % (format(len(seq), ","),
                                                                          format(RESOLVE_TIME_BUDGET_MS, ",")),
                     resolve_request(seq, labels=[lab]),
                     [OK, INV, ex("resolve_runs", label=lab, full=True, unless_stopped=True)],
                     truth=dict(base), tags=["budget", "machine-dependent", "warm-up"], timing_class="slow")
        full["warms"] = sid
        cont_from = len(seq) // 2
        cont = _twin(rb, "R-twin-%04d" % (3 + 2 * i), "the same %s bp from base %s: the continuation" % (
                         format(len(seq), ","), format(cont_from, ",")),
                     resolve_request(seq[cont_from:], labels=[lab]),
                     [OK, INV, ex("resolve_runs", label=lab, full=True, unless_stopped=True)], truth=base,
                     tags=["budget", "machine-dependent"], timing_class="slow")
        out += [full, cont]
        relation(ctx, "prefix_equals", rb["id"], full["id"], "the stopped resolve equals the full resolve on its prefix",
                 labels=[lab], a_from_kmer=0, b_from_kmer=0)
        relation(ctx, "prefix_equals", cont["id"], full["id"],
                 "the resolve from base %d equals the full resolve on those k-mers" % cont_from,
                 labels=[lab], a_from_kmer=cont_from, b_from_kmer=0)
    # zero_equals_absent: an exact count of a forward pattern of L >= k is 0 exactly when /resolve finds a k-mer of
    # it missing from the graph
    eco = ctx.rec("eco")
    pairs = []
    q = st.random_dna(ctx.rng("REL", "zero31"), 31)
    pairs.append(("a random 31-mer", q, {"random": q}, "contexts", [ex("absent")], [ex("resolve_absent")], eco.column))
    q = st.random_dna(ctx.rng("REL", "zero60"), 60)
    pairs.append(("a random 60 bp", q, {"random": q}, "paths", [ex("absent", count="paths")], [ex("resolve_absent")],
                  eco.column))
    src = _by_id(scenarios, "P-exact-0020")
    t = src["truth"]
    pairs.append(("a 31-mer of %s" % t["genome"], t["window"], t, "contexts",
                  [ex("contains_source_context", strand="+", require="all")],
                  [ex("resolve_runs", label=t["taxid"], full=True)], t["taxid"]))
    src = _by_id(scenarios, "L-paths-0007")
    t = src["truth"]
    pairs.append(("60 bp of %s" % t["genome"], t["window"], t, "paths", [ex("path_spells", strand="+")],
                  [ex("resolve_runs", label=t["taxid"], full=True)], t["taxid"]))
    src = _by_id(scenarios, "P-exact-0022")
    t = dict(src["truth"])
    r2 = ctx.rng("REL", "zero-variant")
    variant, edits = st.substitute_bases(t["window"], [15], r2)
    t["edits"], t["variant"] = edits, variant
    pairs.append(("a 31-mer of %s with one substitution at 15 (present or not: the routes agree)" % t["genome"],
                  variant, t, "contexts", [], [], t["taxid"]))
    nid = _ids("N-zero")
    out_zero: List[dict] = []
    for title, q, truth, count, pexp, rexp, lab in pairs:
        long = len(q) > K
        preq = pattern_request([_p(q, "z")], strands="forward", long_search="paths" if long else None)
        tags = ["relation", "negative" if "random" in truth else "positive"] + (["smoke"] if not out_zero else [])
        a = scenario(nid(), "N", "%s searched forward: its count" % title, "POST /pattern", preq, truth,
                     [OK, INV, CC] + pexp, tags=tags)
        b = scenario(nid(), "N", "%s resolved: its k-mers in the graph" % title, "POST /resolve",
                     resolve_request(q, labels=[lab]), truth, [OK, INV] + rexp, tags=tags)
        out_zero.append(a)
        out += [a, b]
        relation(ctx, "zero_equals_absent", a["id"], b["id"], "%s: exact 0 on /pattern iff a k-mer is absent on /resolve"
                 % title, pattern=0, count=count)
    # subset: suffix scope within any_offset (the records' first 20 bases)
    for key in ("lam", "t4"):
        any_ = _by_title(scenarios, "P", "the first 20 bases of %s in any_offset scope" % key)
        suf = _by_title(scenarios, "P", "the first 20 bases of %s in suffix scope" % key)
        if key == "lam":
            _add_tag(any_, "smoke")
            _add_tag(suf, "smoke")
        relation(ctx, "subset", suf["id"], any_["id"], "suffix scope's contexts are among any_offset's (%s)" % key,
                 a_pattern=0, b_pattern=0, key="context")
    # subset: an IUPAC pattern's contexts hold the exact window's
    for sid in ("P-iupac-0001", "P-iupac-0003", "P-iupac-0005"):
        iu = _by_id(scenarios, sid)
        t = iu["truth"]
        tw = _twin(iu, "P-twin-%04d" % (1 + len([x for x in out if x["id"].startswith("P-twin")])),
                   "the exact window of %s" % sid, pattern_request([_p(t["window"])]),
                   [OK, INV, CC, ex("contains_source_context", strand="+", require="all")])
        out.append(tw)
        relation(ctx, "subset", tw["id"], iu["id"], "the exact window's contexts are among its IUPAC pattern's",
                 a_pattern=0, b_pattern=0, key="context")
    # subset: a partial list cut at 2 within the complete list of the same pattern
    relation(ctx, "subset", _by_title(scenarios, "P", "partial with max_contexts 2: cut max_contexts")["id"],
             _by_title(scenarios, "P", "max_contexts above the cap: lowered and stated")["id"],
             "the partial list cut at 2 is within the complete list", a_pattern=0, b_pattern=0, key="context")
    # same_answer on chosen fields: the three capabilities documents agree
    cid = {x["route"]: x["id"] for x in scenarios if x["area"] == "C" and x["route"].startswith("GET /")
           and x["route"].endswith("capabilities")}
    capability_relations(ctx, cid, caps)
    return out


def capability_relations(ctx: Ctx, cid: Dict[str, str], caps: Dict[str, dict]) -> List[dict]:
    """same_answer relations among the three capabilities documents (DESIGN §3 C: release, index identity and
    server instance consistent; the gate block a projection of /pattern/capabilities). |cid| maps a route to its
    scenario id; the gate block's fields are read from the capabilities snapshot."""
    out = []
    out.append(relation(ctx, "same_answer", cid["GET /capabilities"], cid["GET /traverse/capabilities"],
                        "/capabilities and /traverse/capabilities name the same server",
                        paths=[["release", "release"], ["server_instance", "attempts.server_instance"],
                               ["attempts.server_instance", "attempts.server_instance"],
                               ["algorithm_version", "algorithm_version"],
                               ["feature_level", "feature_level"], ["resolve", "resolve"], ["in_ram", "in_ram"],
                               ["content_encodings", "content_encodings"], ["compression_level", "compression_level"],
                               ["deadline_check.chunk_target_ms", "deadline_check.chunk_target_ms"],
                               ["deadline_check.max_uninterruptible_ms", "deadline_check.max_uninterruptible_ms"],
                               ["deadline_check.poll_stride", "deadline_check.poll_stride"]]))
    out.append(relation(ctx, "same_answer", cid["GET /pattern/capabilities"], cid["GET /capabilities"],
                        "/capabilities' pattern block is the /pattern/capabilities document", paths=[["", "pattern"]]))
    gate = caps["traverse_capabilities"]["pattern"]
    full = caps["pattern_capabilities"]
    paths = [[k, "pattern." + k] for k in sorted(gate) if k in full and k != "caps"]
    paths += [["caps." + k, "pattern.caps." + k] for k in sorted(gate.get("caps", {}))]
    out.append(relation(ctx, "same_answer", cid["GET /pattern/capabilities"], cid["GET /traverse/capabilities"],
                        "the gate block of /traverse/capabilities projects /pattern/capabilities", paths=paths))
    return out


# ================================================================================================= v6: heavy /pattern
# The owner's decisions of 2026-10-16 for gating POST /pattern (manifest v6; README "Manifest v6", round 2).
#
# Warm-ups: the /pattern scenarios whose MCP cold call (repeat 0, which the paired loop runs first) timed out (the
# host's 503 deadline relayed by the service) or was time-limited in the accepted v5 baseline run
# 20261015T161813Z_refseq97-81c620ef, all 25 of them, as the owner's decision 3 has it: each gets a warm-up
# (pattern_warm_up), placed right before it (build_manifest.place_warm_ups), so that its MCP arm's first call meets
# a warm host; the v5 figures of each are stated in the warm-up's truth. In 23 of them a later execution (the warm MCP
# repeat or a direct arm) completed; M-protein-0002 and P-heavy-0001 (P-warm-0024 / -0025, added by round 2's
# verification, L4) completed none, warm runs included, and are budget_stress themselves: their warm-ups run and are
# reported outside the gate, as they are. The ids are fixed per target (a warm-up added later takes the next free
# number; no id moves).
PATTERN_WARM_UPS = {
    "M-planted-0007": ("P-warm-0001", "MCP r0 (cold) a relayed 503 deadline at 63.4 s; MCP r1 complete in 3.6 s; direct-local r0 complete in 24.3 s"),
    "M-planted-0008": ("P-warm-0002", "MCP r0 (cold) time-limited at 63.3 s; MCP r1 complete in 3.7 s; direct-local r0 complete in 16.1 s"),
    "M-planted-0009": ("P-warm-0003", "MCP r0 (cold) time-limited at 63.2 s; MCP r1 complete in 3.8 s; direct-local r0 complete in 2.7 s"),
    "M-planted-0013": ("P-warm-0004", "MCP r0 (cold) time-limited at 63.2 s; MCP r1 complete in 3.6 s; direct-local r0 complete in 15.3 s"),
    "M-planted-0014": ("P-warm-0005", "MCP r0 (cold) time-limited at 63.2 s; MCP r1 complete in 3.7 s; direct-local r0 complete in 20.7 s"),
    "M-planted-0015": ("P-warm-0006", "MCP r0 (cold) a relayed 503 deadline at 63.2 s; MCP r1 complete in 3.8 s; direct-local r0 complete in 10.5 s"),
    "M-planted-0016": ("P-warm-0007", "MCP r0 (cold) time-limited at 63.2 s; MCP r1 complete in 4.1 s; direct-local r0 complete in 10.0 s"),
    "M-planted-0017": ("P-warm-0008", "MCP r0 (cold) time-limited at 63.2 s; MCP r1 complete in 3.6 s; direct-local r0 complete in 19.8 s"),
    "M-planted-0018": ("P-warm-0009", "MCP r0 (cold) time-limited at 63.2 s; MCP r1 complete in 3.8 s; direct-local r0 complete in 1.8 s"),
    "M-protein-0001": ("P-warm-0010", "MCP r0 (cold) a relayed 503 deadline at 63.2 s; MCP r1 complete in 19.1 s; direct-local r0 time-limited at 59.8 s"),
    "M-protein-0003": ("P-warm-0011", "MCP r0 (cold) time-limited at 63.3 s; MCP r1 complete in 61.3 s; direct-local r0 a 503 deadline at 60.2 s"),
    "M-protein-0004": ("P-warm-0012", "MCP r0 (cold) a relayed 503 deadline at 63.1 s; MCP r1 complete in 3.7 s; direct-local r0 complete in 32.0 s"),
    "M-protein-0005": ("P-warm-0013", "MCP r0 (cold) time-limited at 63.2 s; MCP r1 complete in 8.4 s; direct-local r0 complete in 17.7 s"),
    "M-protein-0006": ("P-warm-0014", "MCP r0 (cold) a relayed 503 deadline at 63.2 s; MCP r1 complete in 9.1 s; direct-local r0 time-limited at 59.8 s"),
    "M-protein-0007": ("P-warm-0015", "MCP r0 (cold) time-limited at 63.1 s; MCP r1 complete in 6.7 s; direct-local r0 time-limited at 59.8 s"),
    "N-nearmiss-0006": ("P-warm-0016", "MCP r0 (cold) a relayed 503 deadline at 63.1 s; MCP r1 complete in 3.5 s; direct-local r0 complete in 13.3 s"),
    "N-nearmiss-0013": ("P-warm-0017", "MCP r0 (cold) a relayed 503 deadline at 63.2 s; MCP r1 complete in 3.5 s; direct-local r0 complete in 21.6 s"),
    "N-nearmiss-0020": ("P-warm-0018", "MCP r0 (cold) a relayed 503 deadline at 63.1 s; MCP r1 complete in 3.6 s; direct-local r0 complete in 28.3 s"),
    "N-nearmiss-0027": ("P-warm-0019", "MCP r0 (cold) a relayed 503 deadline at 63.1 s; MCP r1 complete in 3.6 s; direct-local r0 complete in 15.1 s"),
    "N-nearmiss-0029": ("P-warm-0020", "MCP r0 (cold) a relayed 503 deadline at 63.2 s; MCP r1 complete in 3.8 s; direct-local r0 complete in 45.2 s"),
    "N-nearmiss-0030": ("P-warm-0021", "MCP r0 (cold) time-limited at 63.2 s; MCP r1 complete in 3.8 s; direct-local r0 complete in 12.1 s"),
    "N-nearmiss-0032": ("P-warm-0022", "MCP r0 (cold) time-limited at 63.4 s; MCP r1 complete in 3.7 s; direct-local r0 complete in 22.4 s"),
    "N-usecase-0008": ("P-warm-0023", "MCP r0 (cold) a relayed 503 deadline at 63.0 s; MCP r1 complete in 3.4 s; direct-local r0 complete in 45.4 s"),
    "M-protein-0002": ("P-warm-0024", "MCP r0 (cold) time-limited at 63.3 s; MCP r1 time-limited at 63.4 s; direct-local r0 time-limited at 59.8 s; no execution completed"),
    "P-heavy-0001": ("P-warm-0025", "MCP r0 (cold) a relayed 503 deadline at 63.0 s; MCP r1 a relayed 503 deadline at 63.1 s; direct-local r0 time-limited at 59.9 s; no execution completed"),
}
PATTERN_WARM_BASELINE = "20261015T161813Z_refseq97-81c620ef (manifest v5, 60,000 ms budgets)"
PATTERN_WARM_STRESS = ("a warm-up's first run meets the cold host by design (it reads the rows its target reads, so "
                       "that the target's MCP call meets them warm): whether it completes within the budget depends "
                       "on the host's cache and load, not on the deployment")


def pattern_warm_up(target: dict, sid: str, facts: str) -> dict:
    """The warm-up |sid| of the /pattern scenario |target| (PATTERN_WARM_UPS): the target's request as it is, run
    right before it (build_manifest.place_warm_ups) on direct-local only (timing_class long, "> 100 s: direct-local
    only": its cold run may take most of the 110,000 ms budget, and neither the MCP arm, whose first call it warms
    the host for, nor direct-remote behind the ingress runs it), with the generic checks only (status, invariants,
    counts_consistent, budget_respected: the property is the target's), tagged budget_stress (PATTERN_WARM_STRESS:
    outside the gate, still run and reported). |facts|: what the target did in the v5 baseline."""
    t = {"warm_up": "the host's rows of %s's search, read right before it, so that its MCP arm's first call (repeat "
                    "0, cold) meets them warm; v5 baseline %s: %s" % (target["id"], PATTERN_WARM_BASELINE, facts),
         "budget_stress": {"reason": PATTERN_WARM_STRESS}}
    expect = [OK, INV] + ([CC] if any(e["check"] == "counts_consistent" for e in target["expect"]) else [])
    return scenario(sid, "P", "warm-up of %s: its request, run right before it on direct-local" % target["id"],
                    "POST /pattern", copy.deepcopy(target["request"]), t, [dict(e) for e in expect],
                    timing_class="long", tags=("warm-up", "budget_stress", "machine-dependent"), warms=target["id"])


def pattern_warm_ups(scenarios: List[dict]) -> List[dict]:
    """The warm-ups of PATTERN_WARM_UPS whose target |scenarios| holds, in the table's order."""
    by = dict((s["id"], s) for s in scenarios)
    return [pattern_warm_up(by[t], sid, facts) for t, (sid, facts) in PATTERN_WARM_UPS.items() if t in by]


# budget_stress (the owner, 2026-10-16): the heavy /pattern scenarios that no budget the ingress admits (110,000 ms)
# makes decisive on a cold host, by the v5 baseline's figures (the readiness study's projection of each stopped
# answer's step rate); outside the gate (gate_policy's budget_stress rule), still run and reported, as v5's long walks.
PATTERN_BUDGET_STRESS = {
    "M-protein-0001": "a decapeptide at max_mismatches 1 (7.8 M steps): in the v5 baseline its cold MCP call ended in "
                      "a relayed 503 deadline and its direct-local cold run was time-limited at 60 s with 17 % of the "
                      "steps done (projected 345 s); only the warm MCP repeat completed (19 s): no budget below the "
                      "ingress's ~120 s cut makes it decisive on a cold host",
    "M-protein-0002": "a decapeptide at max_mismatches 1: not one of its four executions in the v5 baseline finished "
                      "within 60 s, the warm MCP repeat included (time-limited at 59.9 s): its completion depends on "
                      "the host's cache and load at any budget the ingress admits",
    "M-protein-0003": "a decapeptide at max_mismatches 1 (7.8 M steps): in the v5 baseline its cold MCP call did 5 % "
                      "of the steps in 60 s (projected 1,180 s), its direct-local cold run ended in a 503 deadline, "
                      "and even the warm MCP repeat took 56.7 s",
    "M-protein-0007": "a decapeptide at max_mismatches 1 (3.1 M steps): in the v5 baseline its cold MCP call did 19 % "
                      "of the steps in 60 s (projected 317 s) and its direct-local cold run 47 % (projected 126 s); "
                      "only warm runs completed",
    "P-heavy-0001": "a 12-mer in any_offset scope, count (heavy): not decisive at any budget in the v5 baseline: its "
                    "one 200 was a time-limited at_least at 59.9 s, and 5 of its 6 executions were 503 deadlines at "
                    "60.17 to 60.52 s, the work overrunning the 250 ms finalisation reserve (SPEC-pattern-search §6 / "
                    "§7.6 allow that on any host)",
    "N-nearmiss-0030": "a near-miss decapeptide at max_mismatches 1, mode partial (1.3 M steps): in the v5 baseline "
                       "its cold MCP call did 49 % of the steps in 60 s (projected 122 s, past a 110 s budget); its "
                       "direct and warm runs, after that call, completed in 1 to 12 s",
}


def tag_pattern_budget_stress(scenarios: List[dict]) -> List[str]:
    """Tags the scenarios of PATTERN_BUDGET_STRESS budget_stress with their reason (truth.budget_stress.reason, which
    gate_policy's exclusion states); returns the ids tagged."""
    out = []
    for s in scenarios:
        why = PATTERN_BUDGET_STRESS.get(s["id"])
        if why is None:
            continue
        s["truth"] = dict(s.get("truth") or {}, budget_stress={"reason": why})
        s["tags"] = sorted(set(s["tags"]) | {"budget_stress"})
        out.append(s["id"])
    return out


# ================================================================================================= v6: /pattern via MCP
# What the MCP arm's /pattern answer cannot carry (the service's pattern_search status, app/pattern/merge.py, and the
# async session's normaliser, scripts/bench/suite_mcp_normalise.py, both at 9d76035e, read only: the normalised body
# is mode, patterns and index; an answered entry keeps the merged view's id, kind, length, palindromic,
# information_bits, scope, mode, counts, statements, determinism and notes and the database's work, timing, predicate,
# selection, absence_filter and load_ms; verified against the recorded MCP answers of the v5 baseline): a property
# check that reads one of these is marked not implemented yet via MCP at property level (mcp.not_implemented_checks,
# the owner's rule of 2026-10-15), and build_manifest.scenario_level_marks makes it the scenario's mark where it would
# cover every property check. Not marked, because the MCP's answer DOES carry them, per entry: the host's top-level
# predicate block (by_database[].predicate, verbatim on every entry) and timing.load_ms (by_database[].load_ms): the
# runner's normaliser leaves them on the entry where the host states them at the top (MCP_PATTERN_PER_ENTRY), and the
# checks read them there on the MCP arm (checks.MCP_PATTERN_ENTRY_FIELDS: Q-pred-0011 / -0012, P-edge-0009).
MCP_PATTERN_BODY_KEYS = ("mode", "patterns", "index")
MCP_PATTERN_PER_ENTRY = ("predicate", "timing.load_ms")
MCP_PATTERN_ENTRY_NOT_EXPOSED = {
    "pattern": "the entry's pattern text",
    "strands": "the entry's strands",
    "absence_scope": "the entry's absence_scope (the service states the request's scope instead)",
    "anchor_information_bits": "the entry's anchor_information_bits",
    "min_anchor_information_bits": "the entry's min_anchor_information_bits",
    "limits": "the entry's limits",
    "max_mismatches": "an answered entry's max_mismatches (only a refused slot carries it)",
    "mismatch_information_bits": "an answered entry's mismatch_information_bits (only a refused slot carries it)",
    "min_anchor_mismatch_information_bits": "an answered entry's min_anchor_mismatch_information_bits (only a refused "
                                            "slot carries it)",
}


def _expected_paths(args: dict) -> List[str]:
    """The paths an entry_fields / answer_fields expectation reads a value from (equals, present); `absent` reads
    an absence, which an unexposed field shows as well."""
    return list((args.get("equals") or {}).keys()) + list(args.get("present") or [])


def pattern_mcp_unexposed(sc: dict) -> Dict[str, List[str]]:
    """{check: [what the MCP's /pattern answer cannot carry that it reads]} of a scenario the MCP runs
    (pattern_search): an entry_fields path whose first field is in MCP_PATTERN_ENTRY_NOT_EXPOSED (a refused slot's
    mismatch keys are carried), an answer_fields path outside MCP_PATTERN_BODY_KEYS that the entries do not carry
    either (MCP_PATTERN_PER_ENTRY)."""
    if sc.get("route") != "POST /pattern" or (sc.get("mcp") or {}).get("tool") != "pattern_search":
        return {}
    refused = set((e.get("args") or {}).get("pattern", 0) for e in sc["expect"] if e["check"] == "pattern_slot_error")
    out: Dict[str, List[str]] = {}
    for e in sc["expect"]:
        a = e.get("args") or {}
        if e["check"] == "entry_fields":
            slot = a.get("pattern", 0) in refused
            for p in _expected_paths(a):
                head = p.split(".", 1)[0]
                if head in MCP_PATTERN_ENTRY_NOT_EXPOSED and not (slot and "mismatch" in head):
                    out.setdefault("entry_fields", [])
                    if MCP_PATTERN_ENTRY_NOT_EXPOSED[head] not in out["entry_fields"]:
                        out["entry_fields"].append(MCP_PATTERN_ENTRY_NOT_EXPOSED[head])
        elif e["check"] == "answer_fields" and not a.get("when_http"):
            for p in _expected_paths(a):
                if p.split(".", 1)[0] not in MCP_PATTERN_BODY_KEYS and \
                        not any(p == q or p.startswith(q + ".") for q in MCP_PATTERN_PER_ENTRY):
                    what = "the answer's %s (the MCP's body is mode, patterns and index)" % p
                    if what not in out.setdefault("answer_fields", []):
                        out["answer_fields"].append(what)
    return out


# v7 (the service owner's statement, 2026-10-17): the MCP's /pattern answer leaves out the entry's pattern text BY
# DESIGN, not for want of a feature. A check that reads nothing else of what the answer leaves out is classed
# service-side by design at property level (mcp.by_design_checks [{check, case}], lane V5's manifest-side mechanism:
# the mark sets a not_exposed reading of the check aside on the MCP arm, the check stays required on the direct arms)
# instead of not implemented yet via MCP, under the by-design table's case that the scenario names (mcp.by_design)
# for the service's handling of that field: pattern.lowercase for the pattern text (the service checks the text in
# either case and upper-cases it before it is stored and sent, so the host never sees the caller's text from the MCP
# arm; P-edge-0004 reads the host's upper-cased echo). T-edge-0008 is modelled through the table instead: its case
# traverse.duplicate_seeds states an explains_checks entry {check, detail_matches} for the not_exposed answer_fields,
# which the comparison verifies as a by-design pair; the same here would add such an entry to pattern.lowercase, a
# change of the service's own file (service_by_design.json), which the suite does not make. A scenario that names
# none of a field's cases keeps the not-implemented mark (the builder refuses a by_design_checks case outside the
# table).
MCP_PATTERN_ENTRY_BY_DESIGN = {"pattern": ("pattern.lowercase",)}


def pattern_by_design_case(sc: dict, check: str, what: List[str]) -> Optional[str]:
    """The by-design table's case under which the MCP's /pattern answer leaves out everything |check| of |sc| reads
    of what it cannot carry (|what|: pattern_mcp_unexposed's descriptions), by MCP_PATTERN_ENTRY_BY_DESIGN: an
    entry_fields check whose unexposed fields are all by design, each with a case the scenario names; else None."""
    if check != "entry_fields" or not what:
        return None
    named = (sc.get("mcp") or {}).get("by_design") or []
    heads = dict((MCP_PATTERN_ENTRY_NOT_EXPOSED[h], h) for h in MCP_PATTERN_ENTRY_NOT_EXPOSED)
    cases = []
    for w in what:
        case = next((c for c in MCP_PATTERN_ENTRY_BY_DESIGN.get(heads.get(w), ()) if c in named), None)
        if case is None:
            return None
        cases.append(case)
    return cases[0] if len(set(cases)) == 1 else None


def mark_pattern_mcp_fields(scenarios: List[dict], by_design: bool = False) -> Dict[str, List[str]]:
    """Marks every property check of a /pattern scenario that reads what the MCP's answer cannot carry
    (pattern_mcp_unexposed) not implemented yet via MCP (mcp_not_implemented_check, with the fields as its reason),
    or, with |by_design| (manifest v7), service-side by design where the service leaves out by design all it reads
    (pattern_by_design_case: mcp_by_design_check with the case); returns {id: [checks marked]}."""
    out = {}
    for s in scenarios:
        for check, what in pattern_mcp_unexposed(s).items():
            if any(e.get("check") == check for key in ("not_implemented_checks", "by_design_checks")
                   for e in (s["mcp"].get(key) or [])):
                continue
            case = pattern_by_design_case(s, check, what) if by_design else None
            if case:
                mcp_by_design_check(s, check, case)
            else:
                mcp_not_implemented_check(s, check, "the MCP's /pattern answer does not carry %s" % "; ".join(what))
            out.setdefault(s["id"], []).append(check)
    return out
