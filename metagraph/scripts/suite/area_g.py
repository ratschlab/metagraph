"""Area G of manifest v8 (DESIGN-v8 §3-§5, FORMAT-v8 §1-§2 and v2-v2.3): the graphlet processing tools, run five
ways on a corpus of graphlets of refseq33m. areas.gen_G is this module's gen_G.

A scenario of area G is a `POST /traverse` request with `output.detail: graphlet` and one seed, whose strategy is the
EXACT merged strategy the service sends for its MCP mapping (traverse_fetch with a preset and overrides: the preset's
strategy, the overrides deep-merged into it, output {detail: graphlet, timing: true}; DESIGN-v8 §3, the v7 echo of
T-graphlet-0001 in traverse_status.overrides), plus a `graphlet` block (FORMAT §1):
  twin             {strategy}: the request's strategy with output.detail "full" and without bounds.max_memory_mb
                   (TWIN_WITHOUT: F-D4; direct-local only, W5's input);
  preconditions    structural facts frozen from the rehearsal (graphlet_precheck.PRECONDITION_FACTS; a failed seed
                   states {"no_graphlet": true}, FORMAT v2 D3); empty before the rehearsal;
  library_commit   the library commit the manifest is frozen against (a build parameter; FORMAT v2 D1);
  ops              the op plan (op_plan below): every argument named, run-time selectors for every value read from
                   an earlier answer (owner decision D4: {ref, path}, graphlet_select's named selectors, cursor_of);
  ways             {"default": [...]}: the five ways, or W1/W2/W5 for a strategy the MCP cannot express;
  not_exposed      {tool | tool.argument: the owner's reason}: why W3 does not run a library-only tool or argument;
  facts            the service-check facts the suite's reader cannot derive (graphlet_service_checks), when known;
  slot / candidate the corpus slot (INVENTORY §3.2) and the candidate the rehearsal chose for it;
  reference        report-only: the rehearsal's body digest and engine times (DESIGN-v8 §4, D-R4: never the truth).

The corpus (DESIGN-v8 §4, owner decision D1: 12 graphlets, 5 compare pairs): one graphlet per SLOT; per slot several
CANDIDATES, built from the seeds of v7's T scenarios whose v7 answers show the shape the slot needs (the evidence is
quoted per candidate: run 20261019T111815Z_refseq97-1009f77c, the MCP's cold engine time and the direct-local shape)
or drawn here from the truth records. `build_manifest.py --rehearse G` runs every candidate on staging and writes the
corpus file (rehearse_g.py): the candidate chosen per slot, its measured facts (preconditions), the MCP's strategy
echo and the per-op timeouts. Without a corpus file the build is provisional: the first candidate of every slot the
MCP can express, preconditions empty (to be frozen from the rehearsal). Beside the 5 pairs, the rows pair (ROWS_PAIR,
VERIFY-D-R4 R4-1): two rehearsed candidates of one seed whose compare has rows, with a compare W3 pages.

Pure stdlib; the builder never imports the traversal library (the op plans are checked against it in the tests).
"""

from __future__ import annotations

import copy
import hashlib
import json
import re
from typing import Any, Dict, List, Optional, Tuple

import suite_common as sc
from manifest_lib import (TRACE_PRESET, TRIE_PRESET, WALK_PRESET, deep_merge, mcp_traverse, scenario,
                          traverse_request)

AREA = "G"
MAXB = 262144                      # the service's max_bytes ceiling: every W3 whole-list call is designed to fit it
CUT_PAGE = 2048                    # a paging pass's page (the tools' default ceiling; byte-bound pages)
# pages small enough to page on the corpus (VERIFY-D VD8; measured on the rehearsals' documents, W1): a limitations
# listing takes 2-3 pages at 1,024 bytes on every rich graphlet (one row needs up to ~500 bytes: 384 is too small)
LIMS_PAGE = 1024
COMPARE_PAGE = 1024                # cmp-paged: a walks compare with rows (the constrain / annotate pair: 28 rows)
ROW_PAGE = 4096                    # a paging pass over rows with label lists (walk segments, claims, compare)
SERVICE_CAP_CASE = "graphlet.max_bytes_above_service_cap"   # async's case (service_by_design_g.json, 2026-10-20)
ALL_WAYS = ["W1", "W2", "W3", "W4", "W5"]
DIRECT_ONLY_WAYS = ["W1", "W2", "W5"]
LIBRARY_WAYS = ["W1", "W2", "W4", "W5"]
SMOKE_WAYS = ["W1", "W2", "W4"]
# a walks paging pass: rows with a 60 bp tail, pages of 1,200 bytes (several rows each, more than one page on a
# graphlet of 3 walks or more: the page cut and the cursor reuse need a next page)
PAGED_WALKS = {"rank": "support", "spell": "tail", "tail_bp": 60, "max_bytes": 1200}
SUMMARY_CUT_BYTES = 1600           # sum-cut's ceiling before the rehearsal measures the uncut summary (rehearse_g)
SAMPLE_WALKS = 12                  # per-walk sampling on the richest graphlets: the first 12 walks by id per arm
RICHEST = 3                        # DESIGN-v8 §4: per-walk sampling and paging passes on the 3 richest graphlets
HOST_MAX_TIME_MS = 60000           # the host's max_time_ms (snapshot): the effective budget is min(request, it)
HEADROOM_FRACTION = 0.25           # DESIGN-v8 §3: engine time <= 25 % of the effective time budget
LANE_S = {"interactive": 20.0, "heavy": 120.0}   # DESIGN-v8 §3: the service's lanes; every op <= 1/3 of its lane
HEAVY_TOOLS = ("graphlet_compare", "graphlet_export", "graphlet_subtrie", "graphlet_save", "graphlet_load",
               "traverse_continue")

# The arguments each exposed tool takes on the service (its MCP schemas, 2026-10-19), handle / access_token / cursor
# / export_id aside: an op the service runs names every one of them (FORMAT §1: the service's defaults differ)
SERVICE_ARGS = {
    "graphlet_summary": ("arm", "detail", "max_bytes"),
    "graphlet_walks": ("arm", "rank", "n", "min_bp", "label", "spell", "tail_bp", "route_consistent", "max_bytes"),
    "graphlet_walk": ("arm", "walk", "max_bytes"),
    "graphlet_labels": ("arm", "rank", "n", "min_direct_bp", "at_bp", "name", "max_bytes"),
    "graphlet_sequence": ("arm", "walk", "segment", "start", "end", "orientation", "with_seed"),
    "graphlet_compare": ("b", "arm", "mode", "max_bytes"),
    "graphlet_export": ("format", "arm", "walks", "what"),
}
# arguments the library takes that the service does not (INVENTORY §0): only on ops W3 does not run
LIBRARY_ONLY_ARGS = {
    "graphlet_walks": ("coordinates", "budget"),
    "graphlet_labels": ("coordinates", "budget"),
    "graphlet_export": ("path", "other", "budget"),
}
# the library-only tools' arguments (GraphletTools, mcp_tools.py), handle / cursor / budget aside
LIBRARY_TOOL_ARGS = {
    "graphlet_support": ("arm", "walk", "step", "max_bytes"),
    "graphlet_splits": ("arm", "n", "min_labels_before", "max_bytes"),
    "graphlet_claims": ("arm", "min_bp", "route_consistent", "labels", "max_bytes"),
    "graphlet_subtrie": ("labels", "arm", "mode"),
    "graphlet_save": ("path",),
    "graphlet_load": ("path",),
    "graphlet_list": ("max_bytes",),
    "graphlet_free": ("handle",),
    "traverse_continue": ("arm", "walk", "overrides", "execute", "allow_unverified_index", "max_bytes",
                          "reset_branches"),
}
# The owner's reasons W3 does not run an op (DESIGN-v8 §1, §2: decisions of 2026-10-19; FORMAT v2 D2)
OWNER_LIBRARY_ONLY = ("library-only (the owner, 2026-10-19): graphlet_support, graphlet_splits and graphlet_claims "
                      "stay unexposed on MCP")
OWNER_NOT_SERVED = ("library-only (the owner, 2026-10-19): graphlet_subtrie, graphlet_save / graphlet_load and "
                    "traverse_continue are not exposed by the service")
OWNER_SMOKE = ("store management, not exposed by the service (the owner, 2026-10-19): graphlet_list / graphlet_free "
               "are a two-call smoke of the library ways (DESIGN-v8 §2, D-R15)")
REASON_COORDINATES = ("library-only argument: the service's graphlet_walks and graphlet_labels take no coordinates "
                      "(INVENTORY §0)")
REASON_BUDGET = "library-only argument: the service's graphlet tools take no budget (INVENTORY §0)"
REASON_EXPORT = ("library-only argument: the service's graphlet_export takes no path or other and offers no "
                 "what=compare (INVENTORY §0)")
# the owner's decision of 2026-10-19: annotate, output.sequences false and trace with coordinates are direct-only
# scenarios (W1, W2, W5): staging's service runs without the advanced strategy path (service_settings.
# advanced_enabled false: traverse_fetch overrides set only bounds, direction, exhaustive,
# labels.max_labels_per_node and labels.seed_label_kind), so no MCP job can retrieve these graphlets
REASON_DIRECT_ONLY = "not expressible via MCP on staging (advanced strategy path off)"
# slots the owner moved to v9 (2026-10-19): two databases, and a body over 8 MiB (unreachable under the presets: the
# rehearsal's trie candidates stopped at the preset's 1e8 work units with 305,290-byte bodies)
V9_SLOTS = {"over_8mib": "v9 (the owner, 2026-10-19): a body over 8 MiB is unreachable under the presets (the "
                         "rehearsal's trie candidates stopped at 1e8 work units, 305,290 bytes)",
            "two_databases": "v9 (the owner, 2026-10-19)"}
NOT_EXPOSED_TABLE = {
    "graphlet_support": OWNER_LIBRARY_ONLY, "graphlet_splits": OWNER_LIBRARY_ONLY,
    "graphlet_claims": OWNER_LIBRARY_ONLY, "graphlet_subtrie": OWNER_NOT_SERVED, "graphlet_save": OWNER_NOT_SERVED,
    "graphlet_load": OWNER_NOT_SERVED, "traverse_continue": OWNER_NOT_SERVED, "graphlet_list": OWNER_SMOKE,
    "graphlet_free": OWNER_SMOKE, "graphlet_walks.coordinates": REASON_COORDINATES,
    "graphlet_labels.coordinates": REASON_COORDINATES, "graphlet_walks.budget": REASON_BUDGET,
    "graphlet_export.path": REASON_EXPORT, "graphlet_export.other": REASON_EXPORT,
}
PRESETS = {"walk": WALK_PRESET, "trie": TRIE_PRESET, "trace": TRACE_PRESET}

# ===================================================================================================== the corpus
# Slots (INVENTORY §3.2, DESIGN-v8 §4, brief item 6) in manifest order: (slot, scenario id stem, title, conditional)
SLOTS = (
    ("two_arms", "G-arms2", "two arms, branching on both", False),
    ("merges", "G-merge", "constrain with merges (displayed-parent cases, merge-entered labels)", False),
    ("one_arm", "G-arm1", "one arm, branching", False),
    ("cut_lists", "G-cut", "cut label lists (a small labels.max_labels_per_node)", False),
    ("resource_stop", "G-stop", "a deterministic resource stop (bounds.max_work_units, never time)", False),
    ("depth_zero", "G-depth0", "complete_to_bp 0 (bounds.max_steps 1: the no_bins shape)", False),
    ("linear", "G-linear", "a single linear walk", False),
    ("failed_seed", "G-fail", "a failed seed (no graphlet: service no_graphlet)", False),
    ("long_walk", "G-long", "a walk with its seed over 16 kbp (sequence paging)", False),
    ("trace", "G-trace", "a trace retrieval (support trace, on_reconverge keep)", False),
    ("coordinates", "G-coords", "a trace retrieval with record coordinates (direct-only)", False),
    ("annotate", "G-annot", "annotate mode", False),
    ("no_bases", "G-nobases", "output.sequences false (no_bases, compare unknown)", False),
    ("over_8mib", "G-big", "a body over 8 MiB (the service's in_service false; a by-design plan)", True),
    ("two_databases", "G-twodb", "one seed on two databases", True),
)
SLOT_NAMES = tuple(s[0] for s in SLOTS)
# the 3 richest graphlets (DESIGN-v8 §4) unless the corpus names others (rehearsal: most segments)
DEFAULT_RICH = ("two_arms", "merges", "one_arm")
# Compare pairs (DESIGN-v8 §4, owner decision D1): (pair, partner id stem, the corpus slot it pairs with, title)
PAIRS = (
    ("radius", "G-pair-radius", "two_arms", "one seed at two radii (prefix_subset)"),
    ("label_cap", "G-pair-cap", "merges", "one seed under the default and a small label cap (qualified: lower bound)"),
    ("two_seeds", "G-pair-seeds", "two_arms", "two seeds of one locus (comparable false)"),
    ("self", None, "two_arms", "a graphlet with itself"),
    # D1's constrain vs annotate, restored now that annotate is in v8 (direct-only, VERIFY-D VD8): the compare ops run
    # on the annotate graphlet (W1, W2, W5) against the constrain partner, whose own plan runs every way
    ("constrain_annotate", "G-pair-annot", "annotate", "one seed under constrain and under annotate (direct only: "
                                                       "annotate is not expressible via MCP; compare rows)"),
    ("two_databases", None, "two_databases", "one seed on two databases (comparable false)"),
)
EVIDENCE_RUN = "20261019T111815Z_refseq97-1009f77c"


def _by_id(scenarios: List[dict], sid: str) -> dict:
    for s in scenarios:
        if s.get("id") == sid:
            return s
    raise KeyError("no scenario %s" % sid)


def _seed_of(sc: dict) -> dict:
    return copy.deepcopy(sc["request"]["seeds"][0])


def _truth_of(sc: dict) -> dict:
    t = copy.deepcopy(sc.get("truth") or {})
    t.pop("walk", None)   # the source walk's annotation (_annotate_walks) is its radius's, not this graphlet's
    return t


def _cand(cid: str, slot: str, seed: dict, preset: str, overrides: dict, truth: dict, source: str,
          evidence: str, expected: dict, extra: Optional[dict] = None, conditional: bool = False) -> dict:
    """A candidate graphlet: its seed, the preset and the overrides the MCP sends (extra: strategy keys the MCP cannot
    express, making it direct-only), the truth it was drawn from, and what v7's answers showed (expected facts)."""
    return {"cid": cid, "slot": slot, "seed": seed, "preset": preset, "overrides": overrides, "extra": extra or {},
            "truth": truth, "source": source, "evidence": evidence, "expected": expected,
            "conditional": conditional}


def _shift_seed(truth: dict, shift: int) -> Optional[dict]:
    """The window |shift| bases right of |truth|'s (same length) within its flank, or None."""
    flank, fs, s, L = truth.get("flank"), truth.get("flank_start"), truth.get("start"), truth.get("length")
    if not (isinstance(flank, str) and isinstance(fs, int) and isinstance(s, int) and isinstance(L, int)):
        return None
    a = s - fs + shift
    if a < 0 or a + L > len(flank):
        return None
    out = dict(truth, start=s + shift, window=flank[a:a + L])
    return out


def candidates(ctx: Any, built: List[dict]) -> List[dict]:
    """Every candidate of every slot, in slot order (SLOTS), several per slot so that the rehearsal can choose. A
    candidate drawn from a v7 scenario copies its seed and truth; the evidence quotes v7's answers to that seed
    (EVIDENCE_RUN: the MCP's cold engine elapsed ms, then the direct-local graph per arm: walks w, segments s,
    merges m, splits sp, labels)."""
    out: List[dict] = []

    def src(sid: str) -> Tuple[dict, dict]:
        sc = _by_id(built, sid)
        return _seed_of(sc), _truth_of(sc)

    # two arms: derived carrier sets (mpn 70 labels), a header + taxid seed, a column seed
    seed, t = src("T-derived-0003")
    out.append(_cand("two_arms.1", "two_arms", seed, "walk", {"bounds": {"max_extension_bp": 1000}}, t,
                     "T-derived-0003", "mpn derived labels: MCP cold 2399 ms; l w3 s5 sp2, r w4 s22 m5 sp8; 70 labels",
                     {"arms": ["left", "right"], "walks": {"left": 3, "right": 4}, "merges": {"left": 0, "right": 5},
                      "labels": 70}))
    seed, t = src("T-quorum-0002")
    out.append(_cand("two_arms.2", "two_arms", seed, "walk", {"bounds": {"max_extension_bp": 200}}, t,
                     "T-quorum-0002", "mpn header + taxid: MCP cold 636 ms; l w4 s7 sp3, r w3 s8 m1 sp3; 2 labels",
                     {"arms": ["left", "right"], "walks": {"left": 4, "right": 3}, "merges": {"left": 0, "right": 1},
                      "labels": 2}))
    seed, t = src("T-walk-0003")
    out.append(_cand("two_arms.3", "two_arms", seed, "walk", {"bounds": {"max_extension_bp": 100}}, t,
                     "T-walk-0003", "cal column 237561: MCP cold 522 ms; l w2 s3 sp1, r w1 s1; 1 label",
                     {"arms": ["left", "right"], "walks": {"left": 2, "right": 1}, "merges": {"left": 0, "right": 0},
                      "labels": 1}))
    # merges: derived carrier sets with merges on both arms
    seed, t = src("T-derived-0001")
    out.append(_cand("merges.1", "merges", seed, "walk", {"bounds": {"max_extension_bp": 1000}}, t,
                     "T-derived-0001",
                     "buc derived labels: MCP cold 2817 ms; l w2 s27 m8 sp9, r w2 s24 m7 sp8; 8 labels",
                     {"arms": ["left", "right"], "walks": {"left": 2, "right": 2}, "merges": {"left": 8, "right": 7},
                      "labels": 8}))
    seed, t = src("T-derived-0002")
    out.append(_cand("merges.2", "merges", seed, "walk", {"bounds": {"max_extension_bp": 1000}}, t,
                     "T-derived-0002",
                     "mge derived labels: MCP cold 2854 ms; l w1 s4 m1 sp1, r w16 s58 m9 sp24; 6 labels",
                     {"arms": ["left", "right"], "walks": {"left": 1, "right": 16}, "merges": {"left": 1, "right": 9},
                      "labels": 6}))
    # one arm
    seed, t = src("T-walk-0020")
    out.append(_cand("one_arm.1", "one_arm", seed, "walk", {"bounds": {"max_extension_bp": 1000},
                                                           "direction": "left"}, t,
                     "T-walk-0020", "mja header, left: MCP cold 821 ms; l w9 s13 sp4; 1 label",
                     {"arms": ["left"], "walks": {"left": 9}, "merges": {"left": 0}, "labels": 1}))
    seed, t = src("T-derived-0002")
    out.append(_cand("one_arm.2", "one_arm", seed, "walk", {"bounds": {"max_extension_bp": 1000},
                                                           "direction": "right"}, t,
                     "T-derived-0002 (right arm only)", "mge derived labels: both arms MCP cold 2854 ms; r w16 s58 m9",
                     {"arms": ["right"], "walks": {"right": 16}, "merges": {"right": 9}, "labels": 6}))
    seed, t = src("T-walk-0017")
    out.append(_cand("one_arm.3", "one_arm", seed, "walk", {"bounds": {"max_extension_bp": 100}, "direction": "left"},
                     t, "T-walk-0017", "hsa chr17 region header, left: MCP cold 3214 ms; l w4 s7 sp3; 1 label",
                     {"arms": ["left"], "walks": {"left": 4}, "merges": {"left": 0}, "labels": 1}))
    # cut label lists: a small explicit cap, and many labels; a qualified compare of a cut-list graphlet carries the
    # service's hint whatever the cap (graphlet_service_checks.check_hint: each side's own cap, not a fixed 64)
    seed, t = src("T-derived-0001")
    out.append(_cand("cut_lists.1", "cut_lists", seed, "walk",
                     {"bounds": {"max_extension_bp": 1000}, "labels": {"max_labels_per_node": 2}}, t,
                     "T-derived-0001", "buc derived labels (8 labels): MCP cold 2817 ms; cap 2 cuts the lists",
                     {"arms": ["left", "right"], "walks": {"left": 2, "right": 2}, "merges": {"left": 8, "right": 7},
                      "labels": 8, "cut_lists": True}))
    seed, t = src("T-labels-0003")
    out.append(_cand("cut_lists.2", "cut_lists", seed, "walk",
                     {"bounds": {"max_extension_bp": 50},
                      "labels": {"max_labels_per_node": 8, "seed_label_kind": "column"}}, t,
                     "T-labels-0003", "eco column-kind derived labels (873): MCP cold 1094 ms after its warm-up "
                     "(R-warm-0005); l w7 s19 m2, r w18 s37 m1; more than 64 labels per node (hint)",
                     {"arms": ["left", "right"], "walks": {"left": 7, "right": 18}, "merges": {"left": 2, "right": 1},
                      "labels": 873, "cut_lists": True}))
    seed, t = src("T-derived-0003")
    out.append(_cand("cut_lists.3", "cut_lists", seed, "walk",
                     {"bounds": {"max_extension_bp": 1000}, "labels": {"max_labels_per_node": 4}}, t,
                     "T-derived-0003", "mpn derived labels (70): MCP cold 2399 ms; cap 4 cuts the lists, more than 64 "
                     "labels on the seed's nodes (hint)",
                     {"arms": ["left", "right"], "walks": {"left": 3, "right": 4}, "merges": {"left": 0, "right": 5},
                      "labels": 70, "cut_lists": True}))
    # a deterministic resource stop: the work budget above the seed phase (T-budget-0003: 148,723 units)
    seed, t = src("T-budget-0003")
    out.append(_cand("resource_stop.1", "resource_stop", seed, "walk",
                     {"bounds": {"max_extension_bp": 1000, "max_work_units": 300000}}, t,
                     "T-budget-0003", "eco header, work 300,000: MCP cold 450 ms; stop at 35 bp on both arms",
                     {"arms": ["left", "right"], "walks": {"left": 1, "right": 1}, "labels": 1, "stop": "work"}))
    seed, t = src("T-derived-0002")
    out.append(_cand("resource_stop.2", "resource_stop", seed, "walk",
                     {"bounds": {"max_extension_bp": 1000, "max_work_units": 2000000}}, t,
                     "T-derived-0002", "mge derived labels (branching), work 2,000,000: a stop inside the branching "
                     "right arm expected (the cost per bp is not known offline)",
                     {"arms": ["left", "right"], "labels": 6, "stop": "work"}))
    # complete_to_bp 0: one step (no_bins)
    seed, t = src("T-quorum-0002")
    out.append(_cand("depth_zero.1", "depth_zero", seed, "walk",
                     {"bounds": {"max_extension_bp": 200, "max_steps": 1}}, t,
                     "T-quorum-0002", "mpn header + taxid under max_steps 1 (the no_bins fixture's cap)",
                     {"arms": ["left", "right"], "labels": 2, "complete_to_bp": 0}))
    seed, t = src("T-walk-0003")
    out.append(_cand("depth_zero.2", "depth_zero", seed, "walk",
                     {"bounds": {"max_extension_bp": 100, "max_steps": 1}}, t,
                     "T-walk-0003", "cal column under max_steps 1",
                     {"arms": ["left", "right"], "labels": 1, "complete_to_bp": 0}))
    # a single linear walk
    for i, (sid, ext, ms) in enumerate((("T-walk-0034", 1000, 866), ("T-walk-0016", 1000, 1335),
                                        ("T-walk-0001", 100, 407))):
        seed, t = src(sid)
        direction = _by_id(built, sid)["request"]["strategy"]["direction"]
        out.append(_cand("linear.%d" % (i + 1), "linear", seed, "walk",
                         {"bounds": {"max_extension_bp": ext}, "direction": direction}, t, sid,
                         "header, %s: MCP cold %d ms; %s w1 s1; 1 label" % (direction, ms, direction[0]),
                         {"arms": [direction], "walks": {direction: 1}, "segments": {direction: 1}, "labels": 1}))
    # a failed seed: the work budget below the seed phase's cost
    seed, t = src("T-budget-0003")
    out.append(_cand("failed_seed.1", "failed_seed", seed, "walk",
                     {"bounds": {"max_extension_bp": 1000, "max_work_units": 10000}}, t,
                     "T-budget-0003", "eco header: its seed phase took 148,723 units (v4); 10,000 cannot hold it",
                     {"failed": True}))
    seed, t = src("T-walk-0001")
    out.append(_cand("failed_seed.2", "failed_seed", seed, "walk",
                     {"bounds": {"max_extension_bp": 100, "max_work_units": 1000}}, t,
                     "T-walk-0001", "mtb header, 31 bp seed: 1,000 units below any seed phase (G18: >= 435,722 for "
                     "a lambda seed)", {"failed": True}))
    # a walk with its seed over 16 kbp: a long seed and a 10 kbp arm on a single-record column (rows only)
    out += _long_walk_candidates(ctx)
    # trace (the service's trace preset: no coordinates)
    for i, sid in enumerate(("T-coords-0002", "T-coords-0003")):
        seed, t = src(sid)
        out.append(_cand("trace.%d" % (i + 1), "trace", seed, "trace", {"bounds": {"max_extension_bp": 100}}, t,
                         sid, "the trace preset at 100 bp (v7's 250 bp trace + coordinates took 1.7-4.6 s warm on "
                              "direct-local; the MCP cannot set output.coordinates)",
                         {"arms": ["left", "right"], "labels": 1}))
    # trace with record coordinates (direct-only: output.coordinates is not expressible via MCP)
    for i, sid in enumerate(("T-coords-0002", "T-coords-0003")):
        seed, t = src(sid)
        out.append(_cand("coordinates.%d" % (i + 1), "coordinates", seed, "trace",
                         {"bounds": {"max_extension_bp": 100}}, t, sid,
                         "trace + record coordinates at 100 bp (the rehearsal of 2026-10-19 ran the first as trace.3: "
                         "every verdict passed on the direct arm)",
                         {"arms": ["left", "right"], "labels": 1, "coordinates": True},
                         extra={"output": {"coordinates": True}}))
    # annotate (direct-only)
    for i, (sid, ev, exp) in enumerate((
            ("T-labels-0001", "buc annotate: direct-local 1220 ms (warm); l w19 s31 sp12, r w13 s25 sp12; 23 labels",
             {"arms": ["left", "right"], "walks": {"left": 19, "right": 13}, "labels": 23}),
            ("T-labels-0002", "mge annotate: direct-local 551 ms (warm); l w1, r w2 s3 sp1; 5 labels",
             {"arms": ["left", "right"], "walks": {"left": 1, "right": 2}, "labels": 5}))):
        seed, t = src(sid)
        out.append(_cand("annotate.%d" % (i + 1), "annotate", seed, "walk", {"bounds": {"max_extension_bp": 200}}, t,
                         sid, ev, dict(exp, annotate=True),
                         extra={"labels": {"mode": "annotate"}, "branching": {"max_label_branches": "unlimited"}}))
    # output.sequences false (direct-only)
    for i, sid in enumerate(("T-derived-0001", "T-derived-0003")):
        seed, t = src(sid)
        out.append(_cand("no_bases.%d" % (i + 1), "no_bases", seed, "walk", {"bounds": {"max_extension_bp": 1000}},
                         t, sid, "the seed of %s without bases (direct-only: output.sequences is not expressible via "
                                 "MCP)" % sid, {"arms": ["left", "right"], "bases": False},
                         extra={"output": {"sequences": False}}))
    # a body over 8 MiB (conditional): a trie of eco's column-kind derived carrier set, lists uncut
    seed, t = src("T-labels-0003")
    for i, ext in enumerate((500, 1000)):
        out.append(_cand("over_8mib.%d" % (i + 1), "over_8mib", seed, "trie",
                         {"bounds": {"max_extension_bp": ext},
                          "labels": {"max_labels_per_node": 10000, "seed_label_kind": "column"}}, t,
                         "T-labels-0003", "eco column-kind derived (873 labels), the trie preset, uncut lists, %d bp: "
                         "a body over 8 MiB hoped for; v7's 50 bp walk read 6-12 bp/s cold in v4 (headroom unlikely: "
                         "else v9, DESIGN-v8 §4)" % ext, {"over_8mib": True}, conditional=True))
    return out


LONG_SEED_BP = 6600
LONG_RADIUS_BP = 10000


def _long_walk_candidates(ctx: Any) -> List[dict]:
    """A 6,600 bp seed of a single-record column (its rows carry few columns: no coordinates to read) with a 10,000 bp
    right arm (the service maximum): seed + walk 16,600 bases > the 16,384-byte sequence ceiling, so
    graphlet_sequence(with_seed) pages (truncated + next_from). Drawn here (own random stream: no other area's
    draws move); the walk preset (15 s) and the trie preset (30 s) for the rehearsal's headroom."""
    out = []
    # .4 / .5 (2026-10-19, after the first rehearsal had walked the windows of .1-.3): fresh windows of the same
    # records, so that the second rehearsal measures the warm-ups' effect from a cold host
    for i, (key, preset, stream) in enumerate((("ser", "walk", ""), ("mpn", "walk", ""), ("ser", "trie", ""),
                                                ("ser", "walk", "fresh"), ("mpn", "walk", "fresh"))):
        try:
            rec = ctx.rec(key)
            rng = ctx.rng("G", "long_walk", key, stream) if stream else ctx.rng("G", "long_walk", key)
            t = _long_window(ctx, rng, rec, LONG_SEED_BP, LONG_RADIUS_BP + 200)
        except (KeyError, ValueError) as e:
            ctx.skip(AREA, "long_walk.%d" % (i + 1), "no %d bp window: %s" % (LONG_SEED_BP, e))
            continue
        seed = {"sequence": t.pop("window"), "labels": [rec.column]}
        t["flank"] = ""
        t["window_sha256"] = hashlib.sha256(seed["sequence"].encode()).hexdigest()
        out.append(_cand("long_walk.%d" % (i + 1), "long_walk", seed, preset,
                         {"bounds": {"max_extension_bp": LONG_RADIUS_BP}, "direction": "right"}, t,
                         "drawn (%s, column %s, %d record)" % (key, rec.column, rec.records_in_column),
                         "a %d bp seed + a %d bp right arm under the %s preset; cold cost unknown (~1 ms/bp cold for "
                         "v7's single-record walks: the rehearsal decides)" % (LONG_SEED_BP, LONG_RADIUS_BP, preset),
                         {"arms": ["right"], "labels": 1, "long": True}))
    return out


def _long_window(ctx: Any, rng: Any, rec: Any, L: int, margin: int, tries: int = 400) -> dict:
    """A window of L bases of |rec| at least |margin| bases from its ends with a right flank of |margin| bases, ACGT
    only, unique on both strands of the record: Ctx.window's rules without its record of windows drawn before, so
    that the builder (after every other area's draws) and the rehearsal (a fresh context) draw the same window."""
    import seqtools as st
    seq = rec.seq
    lo, hi = margin, len(seq) - L - margin
    if rec.region or hi < lo:
        raise ValueError("no room for a %d bp window with %d bp flanks in %s" % (L, margin, rec.key))
    for _ in range(tries):
        s = rng.randint(lo, hi)
        w = seq[s:s + L]
        right = seq[s + L:s + L + margin]
        if not st.is_acgt(w) or not st.is_acgt(right) or st.low_complexity_dust_score(w[:64]) > 2.0:
            continue
        if seq.count(w) != 1 or st.revcomp(w) in seq:
            continue
        return ctx.truth_at(rec, s, L)
    raise ValueError("no unique %d bp window found in %s" % (L, rec.key))


# The full twin's one difference from the request's strategy (the coordinator's decision on F-D4, 2026-10-19):
# SPEC-core §7.1 R21 (4) keeps a merge's arrival order in `full` detail under a memory budget, while `graphlet` detail
# shows its majority parent first, so with the presets' max_memory_mb the twin's first-parent paths differ from the
# graphlet's wherever a merge's majority arrived second (the rehearsal of 2026-10-19: 9 candidates). The twin is
# direct-local only and feeds W5; the budget does not bind at these sizes. graphlet_precheck allows exactly this key.
TWIN_WITHOUT = ("bounds", "max_memory_mb")


def twin_strategy(strategy: dict) -> dict:
    """The full twin's strategy: |strategy| with output.detail "full" and without bounds.max_memory_mb
    (TWIN_WITHOUT)."""
    twin = copy.deepcopy(strategy)
    twin.setdefault("output", {})["detail"] = "full"
    (twin.get(TWIN_WITHOUT[0]) or {}).pop(TWIN_WITHOUT[1], None)
    return twin


# The owner's decision of 2026-10-19: the long walk gets warm-up resolves of its region right before it, as the P-warm
# scenarios do (direct-local only, timing_class long: the MCP arm's first call meets warm rows). The rehearsal of that
# day: its cold MCP walks took 14.7 s and 13.1 s against the walk preset's 15 s (25 % is 3.75 s), warm ~0.2 s.
LONG_WARM_IDS = ("G-warm-0001", "G-warm-0002")
LONG_WARM_MAX_QUERY_BP = 10000     # the host's max_query_bp (its /traverse/capabilities): a direct-local resolve's cap
LONG_WARM_MARGIN_BP = 45           # areas.WARM_MARGIN_BP: what a warm-up reads beyond the walk


def long_walk_warm_ups(ctx: Any, c: dict, target_id: str, ids: Tuple[str, ...] = LONG_WARM_IDS) -> List[dict]:
    """The warm-ups of the long walk |c| (LONG_WARM_IDS): its seed and its right arm's radius plus
    LONG_WARM_MARGIN_BP, in two contiguous pieces of at most LONG_WARM_MAX_QUERY_BP overlapping by k - 1, each resolved
    against the walk's own label (a column: its rows), every k-mer supported; direct-local only (timing_class long, the
    MCP mapping not implemented: bounds.max_query_bp above the service's 2,000), budget_stress and machine-dependent
    as the P-warm warm-ups (their own cold run meets the cold host by design)."""
    from manifest_lib import RESOLVE_TIME_BUDGET_MS, resolve_request
    t = c["truth"]
    rec = ctx.rec(t["genome"])
    label = c["seed"]["labels"][0]
    radius = c["overrides"]["bounds"]["max_extension_bp"]
    s0 = t["start"] - rec.offset
    a, b = s0, min(len(rec.seq), s0 + t["length"] + radius + LONG_WARM_MARGIN_BP)
    k = 31
    n = len(ids)
    piece = -(-(b - a + (n - 1) * (k - 1)) // n)
    assert piece <= LONG_WARM_MAX_QUERY_BP, "a long-walk warm-up within the host's max_query_bp"
    out = []
    for i, sid in enumerate(ids):
        start = a + i * (piece - (k - 1))
        end = min(b, start + piece) if i + 1 < n else b
        region = ctx.truth_at(rec, start, end - start, flank=0)
        region["flank"] = ""
        query = region.pop("window")
        region["window_sha256"] = hashlib.sha256(query.encode()).hexdigest()
        region["warm_up"] = ("the rows of %s's walk (piece %d of %d: its %s bp seed and %s bp right arm + %d), read "
                             "right before it on direct-local so that its MCP arm's first call meets them warm; the "
                             "rehearsal of 2026-10-19: cold MCP walks 14.7 s / 13.1 s against 15 s, warm ~0.2 s" % (
                                 target_id, i + 1, n, format(t["length"], ","), format(radius, ","),
                                 LONG_WARM_MARGIN_BP))
        region["budget_stress"] = {"reason": "a warm-up's first run meets the cold host by design"}
        req = resolve_request(query, labels=[label], time_budget_ms=RESOLVE_TIME_BUDGET_MS,
                              max_query_bp=LONG_WARM_MAX_QUERY_BP)
        out.append(scenario(sid, AREA, "warm-up of %s: %s bp of its region resolved against %s on direct-local" % (
                                target_id, format(len(query), ","), label), "POST /resolve", req, region,
                            [{"check": "status", "args": {"http": 200}}, {"check": "invariants", "args": {}},
                             {"check": "resolve_runs", "args": {"label": label, "full": True}}],
                            timing_class="long", tags=("warm-up", "budget_stress", "machine-dependent", "area-g"),
                            warms=target_id))
    return out


def candidate_strategy(c: dict) -> dict:
    """The exact strategy of a candidate's request: the preset's, the overrides deep-merged, output {detail:
    graphlet, timing: true}, then the direct-only keys (manifest_lib.traverse_request)."""
    return traverse_request([c["seed"]], c["overrides"], preset=c["preset"], detail="graphlet",
                            extra=c["extra"] or None)["strategy"]


def mcp_merge(preset: str, overrides: Optional[dict]) -> dict:
    """What the service sends for traverse_fetch(preset, overrides): the preset's strategy with the overrides
    deep-merged and output {detail: graphlet, timing: true} (the v7 echo of T-graphlet-0001)."""
    return deep_merge(deep_merge(PRESETS[preset], overrides or {}), {"output": {"detail": "graphlet", "timing": True}})


def expressible(c: dict) -> bool:
    """Whether the MCP runs the candidate: a strategy the service's presets and overrides express, and no reason of
    the corpus to keep it direct-only (direct_only_reason: the constrain partner of the direct-only annotate pair)."""
    return not c["extra"] and not c.get("direct_only_reason")


def direct_reason(c: dict) -> str:
    return c.get("direct_only_reason") or REASON_DIRECT_ONLY


REASON_ANNOTATE_PAIR = ("the constrain / annotate pair runs on the direct ways (annotate is not expressible via MCP on "
                        "staging: advanced strategy path off)")


# ===================================================================================================== op plans
def ref(op: str, path: str) -> dict:
    return {"ref": op, "path": path}


def sel(name: str, **args: Any) -> dict:
    return {"select": name, "args": args}


def _op(oid: str, tool: str, args: dict, paging: str = "none", ways: Optional[List[str]] = None, **fields: Any) -> dict:
    o = {"id": oid, "tool": tool, "args": args, "paging": paging}
    if ways is not None:
        o["ways"] = list(ways)
    o.update(fields)
    return o


def W(arm: Optional[str], **kw: Any) -> dict:
    return dict({"arm": arm, "rank": "id", "n": 10000, "min_bp": 0, "label": None, "spell": "none", "tail_bp": 60,
                 "route_consistent": True, "max_bytes": MAXB}, **kw)


def L(**kw: Any) -> dict:
    return dict({"arm": None, "rank": "direct_bp", "n": 10000, "min_direct_bp": 0, "at_bp": None, "name": None,
                 "max_bytes": MAXB}, **kw)


def S(arm: Optional[str], **kw: Any) -> dict:
    return dict({"arm": arm, "walk": None, "segment": None, "start": 0, "end": None, "orientation": "natural",
                 "with_seed": False}, **kw)


def E(**kw: Any) -> dict:
    return dict({"format": "fasta", "arm": None, "walks": None, "what": None}, **kw)


def C(b: str, **kw: Any) -> dict:
    return dict({"b": {"graphlet": b}, "arm": None, "mode": "claims", "max_bytes": MAXB}, **kw)


def lib(o: dict, reason: str) -> dict:
    """|o| as a library-only op: the library ways and W5, W3 not_exposed with the owner's reason."""
    o["ways"] = list(LIBRARY_WAYS)
    o["not_exposed"] = {"W3": reason}
    return o


# A paging pass (an op id holding "-paged") must really page on its graphlet (the coordinator, round 5): its
# max_bytes is the graphlet's measured one, facts page_bytes {op id: bytes}, the largest page at which the local
# library (W1) took at least two pages of the listing on the rehearsal's document (D-scratch/r5/measure_pages.py);
# null: the listing has one row and cannot page there, so the pass is left out; without a measurement the default.
def page_of(facts: dict, oid: str, default: int) -> Optional[int]:
    pb = facts.get("page_bytes")
    if isinstance(pb, dict) and oid in pb:
        v = pb[oid]
        return v if isinstance(v, int) and not isinstance(v, bool) and v >= 64 else None
    return default


class Plan(object):
    """An op plan under construction: unique op ids, the facts of the graphlet."""

    def __init__(self, sid: str, facts: dict, rich: bool):
        self.sid = sid
        self.facts = facts
        self.rich = rich
        self.ops: List[dict] = []
        self.ids = set()

    def add(self, o: dict) -> dict:
        if o["id"] in self.ids:
            raise ValueError("%s: duplicate op id %s" % (self.sid, o["id"]))
        self.ids.add(o["id"])
        self.ops.append(o)
        return o

    def arms(self) -> List[str]:
        return list(self.facts.get("arms") or ["left", "right"])

    def walks(self, arm: str) -> int:
        return int((self.facts.get("walks") or {}).get(arm) or 0)

    def merges(self, arm: str) -> int:
        return int((self.facts.get("merges") or {}).get(arm) or 0)

    def labels(self) -> int:
        return int(self.facts.get("labels") or 0)


def _arm_arg(p: Plan, arm: str) -> Optional[str]:
    return arm


def plan_summary(p: Plan) -> None:
    arms = p.arms()
    p.add(_op("sum", "graphlet_summary", {"arm": None, "detail": "summary", "max_bytes": MAXB}))
    for a in arms:
        p.add(_op("sum-%s" % a, "graphlet_summary", {"arm": a, "detail": "summary", "max_bytes": MAXB}))
    p.add(_op("lims", "graphlet_summary", {"arm": None, "detail": "limitations", "max_bytes": MAXB}, paging="all"))
    lims = page_of(p.facts, "lims-paged", LIMS_PAGE)
    if p.rich and lims:
        p.add(_op("lims-paged", "graphlet_summary", {"arm": None, "detail": "limitations", "max_bytes": lims},
                  paging="all"))
    p.add(_op("sum-bad-detail", "graphlet_summary", {"arm": None, "detail": "full", "max_bytes": 2048}))
    p.add(_op("sum-bad-arm", "graphlet_summary", {"arm": "up", "detail": "summary", "max_bytes": 2048}))
    if len(arms) == 1:
        other = "left" if arms[0] == "right" else "right"
        p.add(_op("sum-arm-missing", "graphlet_summary", {"arm": other, "detail": "summary", "max_bytes": 2048}))


def walk_ops_of(p: Plan, a: str) -> List[str]:
    """The graphlet_walk op ids of arm |a| (for the walk selectors), in order."""
    return [o["id"] for o in p.ops if o["tool"] == "graphlet_walk" and o["args"].get("arm") == a
            and "max_bytes" in o["args"] and o["args"]["max_bytes"] == MAXB and "handle" not in o["args"]
            and not o["id"].endswith(("-out", "-neg"))]


def plan_walks_and_walk(p: Plan) -> None:
    arms = p.arms()
    for a in arms:
        n = p.walks(a)
        p.add(_op("walks-id-%s" % a, "graphlet_walks", W(a), paging="all"))
        for rank in ("support", "length", "loss"):
            p.add(_op("walks-%s-%s" % (rank, a), "graphlet_walks", W(a, rank=rank), paging="all"))
        p.add(_op("walks-tail-%s" % a, "graphlet_walks", W(a, rank="length", spell="tail", tail_bp=60), paging="all"))
        p.add(_op("walks-tail-long-%s" % a, "graphlet_walks", W(a, rank="length", spell="tail", tail_bp=20000),
                  paging="all"))
        p.add(_op("walks-full-%s" % a, "graphlet_walks", W(a, rank="length", spell="full"), paging="all"))
        mid = max(0, n // 2)
        p.add(_op("walks-minbp-mid-%s" % a, "graphlet_walks",
                  W(a, rank="length", min_bp=ref("walks-length-%s" % a, "/rows/%d/length_bp" % mid)), paging="all"))
        p.add(_op("walks-minbp-past-%s" % a, "graphlet_walks",
                  W(a, rank="length", min_bp=sel("slice_past_end", length_of="walks-length-%s" % a, part="end",
                                                 walk=ref("walks-length-%s" % a, "/rows/0/walk"), past=1)),
                  paging="all"))
        # graphlet_walk: the sample (every walk up to SAMPLE_WALKS on the richest, else first / last / longest)
        p.add(_op("walk-longest-%s" % a, "graphlet_walk",
                  {"arm": a, "walk": ref("walks-length-%s" % a, "/rows/0/walk"), "max_bytes": MAXB}))
        idx = list(range(min(n, SAMPLE_WALKS))) if p.rich else sorted(set([0] + ([n - 1] if n > 1 else [])))
        for i in idx:
            p.add(_op("walk-%d-%s" % (i, a), "graphlet_walk",
                      {"arm": a, "walk": ref("walks-id-%s" % a, "/rows/%d/walk" % i), "max_bytes": MAXB}))
        wpage = page_of(p.facts, "walk-paged-%s" % a, ROW_PAGE)
        if p.rich and wpage:
            p.add(_op("walk-paged-%s" % a, "graphlet_walk",
                      {"arm": a, "walk": ref("walks-length-%s" % a, "/rows/0/walk"), "max_bytes": wpage},
                      paging="all"))
        wspage = page_of(p.facts, "walks-paged-%s" % a, PAGED_WALKS["max_bytes"])
        if p.rich and wspage:
            p.add(_op("walks-paged-%s" % a, "graphlet_walks", W(a, **dict(PAGED_WALKS, max_bytes=wspage)),
                      paging="all"))
        # the rehearsal's shadow run: merge_walk False when no sampled walk passes a merge
        mw = (p.facts.get("merge_walk") or {}).get(a)
        if p.merges(a) and (mw if mw is not None else p.rich):
            p.add(_op("walk-merge-%s" % a, "graphlet_walk",
                      {"arm": a, "walk": sel("walk_through_merge", walks=walk_ops_of(p, a), arm=a),
                       "max_bytes": MAXB}))
        if p.facts.get("continuation", {}).get(a):
            p.add(_op("walk-cont-%s" % a, "graphlet_walk",
                      {"arm": a, "walk": sel("walk_with_continuation", walks=walk_ops_of(p, a), arm=a),
                       "max_bytes": MAXB}))
    first = arms[0]
    # errors (arm-independent: on the first arm)
    if len(arms) == 2:
        p.add(_op("walks-arm-none", "graphlet_walks", W(None)))           # both arms present: bad_arm
    else:
        p.add(_op("walks-arm-only", "graphlet_walks", W(None), paging="all"))   # the only arm
    p.add(_op("walks-bad-rank", "graphlet_walks", W(first, rank="best")))
    p.add(_op("walks-n0", "graphlet_walks", W(first, n=0)))
    p.add(_op("walks-unknown-label", "graphlet_walks", W(first, label={"name": "no such label"})))
    p.add(_op("walk-out", "graphlet_walk", {"arm": first, "walk": ref("walks-id-%s" % first, "/total"),
                                            "max_bytes": 2048}))
    p.add(_op("walk-neg", "graphlet_walk", {"arm": first, "walk": -1, "max_bytes": 2048}))


def plan_labels(p: Plan) -> None:
    arms = p.arms()
    for arm in [None] + arms:
        tag = arm or "all"
        for rank in ("direct_bp", "reach_bp"):
            p.add(_op("labels-%s-%s" % (tag, rank.split("_")[0]), "graphlet_labels", L(arm=arm, rank=rank),
                      paging="all"))
    p.add(_op("labels-min", "graphlet_labels", L(min_direct_bp=ref("labels-all-direct", "/rows/0/direct_bp")),
              paging="all"))
    for a in arms:
        p.add(_op("labels-at0-%s" % a, "graphlet_labels", L(arm=a, at_bp=0), paging="all"))
        mid = max(1, int(p.facts.get("radius") or 2) // 2)
        p.add(_op("labels-atmid-%s" % a, "graphlet_labels", L(arm=a, at_bp=mid), paging="all"))
        p.add(_op("labels-atend-%s" % a, "graphlet_labels",
                  L(arm=a, at_bp=ref("walks-length-%s" % a, "/rows/0/length_bp")), paging="all"))
        p.add(_op("labels-atpast-%s" % a, "graphlet_labels",
                  L(arm=a, at_bp=sel("slice_past_end", length_of="walks-length-%s" % a, part="end",
                                     walk=ref("walks-length-%s" % a, "/rows/0/walk"), past=1)), paging="all"))
    lpage = page_of(p.facts, "labels-paged", 1024)
    if p.rich and lpage:
        p.add(_op("labels-paged", "graphlet_labels", L(max_bytes=lpage), paging="all"))
    # name=: up to 3 labels, by {ref}, by {name}, by a bare string; per arm and both
    n = min(3, max(1, p.labels()))
    forms = [("ref", lambda i: sel("top_label", labels="labels-all-direct", index=i)),
             ("name", lambda i: {"name": sel("top_label", labels="labels-all-direct", index=i, form="name")}),
             ("bare", lambda i: sel("top_label", labels="labels-all-direct", index=i, form="ref"))]
    for i in range(n):
        form, mk = forms[i % len(forms)]
        p.add(_op("labels-name%d-%s" % (i, form), "graphlet_labels", L(name=mk(i)), paging="all"))
    for a in arms:
        p.add(_op("labels-name0-%s" % a, "graphlet_labels",
                  L(arm=a, name=sel("top_label", labels="labels-%s-direct" % a)), paging="all"))
    p.add(_op("labels-unknown", "graphlet_labels", L(name={"name": "no such label"})))


def plan_sequence(p: Plan) -> None:
    arms = p.arms()
    for a in arms:
        walks = [("longest", ref("walks-length-%s" % a, "/rows/0/walk"))]
        if p.rich and p.walks(a) > 1:
            walks.append(("first", ref("walks-id-%s" % a, "/rows/0/walk")))
        for wname, w in walks:
            for orient in ("natural", "walk"):
                for seed in (False, True):
                    p.add(_op("seq-%s-%s-%s%s" % (wname, a, orient, "-seed" if seed else ""), "graphlet_sequence",
                              S(a, walk=w, orientation=orient, with_seed=seed), paging="all"))
        for orient in ("natural", "walk"):
            p.add(_op("seq-seg-%s-%s" % (a, orient), "graphlet_sequence",
                      S(a, segment=ref("walk-longest-%s" % a, "/rows/0/segment"), orientation=orient),
                      paging="all"))
    a = arms[0]
    w = ref("walks-length-%s" % a, "/rows/0/walk")
    base = "seq-longest-%s-natural" % a

    def end(part: str, **kw: Any) -> dict:
        return sel("slice_past_end", length_of=base, part=part, **kw)
    p.add(_op("seq-slice-in", "graphlet_sequence", S(a, walk=w, start=end("start", keep=20), end=end("start", keep=5)),
              paging="all"))
    p.add(_op("seq-slice-past", "graphlet_sequence", S(a, walk=w, start=end("start", keep=20), end=end("end", past=50)),
              paging="all"))
    p.add(_op("seq-slice-at-length", "graphlet_sequence", S(a, walk=w, start=end("start", keep=0)), paging="all"))
    # errors
    p.add(_op("seq-both", "graphlet_sequence", S(a, walk=w, segment=0)))
    p.add(_op("seq-neither", "graphlet_sequence", S(a)))
    p.add(_op("seq-end-lt-start", "graphlet_sequence", S(a, walk=w, start=10, end=5)))
    p.add(_op("seq-start-past", "graphlet_sequence", S(a, walk=w, start=end("end", past=1))))
    p.add(_op("seq-seed-seg", "graphlet_sequence", S(a, segment=0, with_seed=True)))
    p.add(_op("seq-seg-oor", "graphlet_sequence", S(a, segment=100000)))


def plan_export(p: Plan) -> None:
    arms = p.arms()
    big = (p.facts.get("body_bytes") or 0) > 512 * 1024   # the service's json / gfa body gate (op_memory)
    p.add(_op("exp-mgt", "graphlet_export", E(format="mgt")))
    p.add(_op("exp-fasta", "graphlet_export", E()))
    for a in arms:
        p.add(_op("exp-fasta-%s" % a, "graphlet_export", E(arm=a)))
    a = arms[0]
    walks = [ref("walks-id-%s" % a, "/rows/0/walk"), ref("walks-length-%s" % a, "/rows/0/walk")]
    p.add(_op("exp-fasta-walks", "graphlet_export", E(arm=a, walks=walks)))
    if not big:
        p.add(_op("exp-gfa", "graphlet_export", E(format="gfa")))
        p.add(_op("exp-json", "graphlet_export", E(format="json")))
        p.add(_op("exp-json-walk", "graphlet_export", E(format="json", what="walk", arm=a,
                                                         walks=[ref("walks-length-%s" % a, "/rows/0/walk")])))
    p.add(_op("exp-bad-what", "graphlet_export", E(what="walk")))
    p.add(lib(_op("exp-path", "graphlet_export", dict(E(), path="../escape.fa", other=None)), REASON_EXPORT))


def plan_library_only(p: Plan, partner: Optional[str]) -> None:
    """support, splits, claims (owner: library-only), coordinates and budget arguments, traverse_continue."""
    arms = p.arms()
    for a in arms:
        wl = ref("walks-length-%s" % a, "/rows/0/walk")
        for step in (None, 50):
            p.add(lib(_op("support-longest-%s-%s" % (a, step or "all"), "graphlet_support",
                          {"arm": a, "walk": wl, "step": step, "max_bytes": MAXB}, paging="all"), OWNER_LIBRARY_ONLY))
        if any(o["id"] == "walk-merge-%s" % a for o in p.ops):
            p.add(lib(_op("support-merge-%s" % a, "graphlet_support",
                          {"arm": a, "walk": ref("walk-merge-%s" % a, "/walk"), "step": None, "max_bytes": MAXB},
                          paging="all"), OWNER_LIBRARY_ONLY))
        for mlb in (0, 2, 100000):
            p.add(lib(_op("splits-%s-%d" % (a, mlb), "graphlet_splits",
                          {"arm": a, "n": 10000, "min_labels_before": mlb, "max_bytes": MAXB}, paging="all"),
                      OWNER_LIBRARY_ONLY))
        for mb_name, mb in (("0", 0), ("mid", ref("labels-%s-direct" % a, "/rows/0/direct_bp"))):
            for rc in (True, False):
                p.add(lib(_op("claims-%s-%s-%s" % (a, mb_name, "rc" if rc else "all"), "graphlet_claims",
                              {"arm": a, "min_bp": mb, "route_consistent": rc, "labels": None, "max_bytes": MAXB},
                              paging="all"), OWNER_LIBRARY_ONLY))
        p.add(lib(_op("claims-%s-one" % a, "graphlet_claims",
                      {"arm": a, "min_bp": 0, "route_consistent": True,
                       "labels": sel("top_label", labels="labels-%s-direct" % a), "max_bytes": MAXB}, paging="all"),
                  OWNER_LIBRARY_ONLY))
        if p.labels() >= 2:
            p.add(lib(_op("claims-%s-list" % a, "graphlet_claims",
                          {"arm": a, "min_bp": 0, "route_consistent": True,
                           "labels": [sel("top_label", labels="labels-all-direct", index=0),
                                      sel("top_label", labels="labels-all-direct", index=1)],
                           "max_bytes": MAXB}, paging="all"), OWNER_LIBRARY_ONLY))
        if any(o["id"] == "walk-cont-%s" % a for o in p.ops):
            # a walk that ended for a semantic reason has no continuation (bad_argument): only a continuation walk
            p.add(lib(_op("continue-cont-%s" % a, "traverse_continue",
                          {"arm": a, "walk": ref("walk-cont-%s" % a, "/walk"), "overrides": None, "execute": False,
                           "allow_unverified_index": False, "max_bytes": 16384, "reset_branches": False}),
                      OWNER_NOT_SERVED))
    a = arms[0]
    spage, cpage = page_of(p.facts, "splits-paged", 1024), page_of(p.facts, "claims-paged", ROW_PAGE)
    if p.rich and spage:
        p.add(lib(_op("splits-paged", "graphlet_splits", {"arm": a, "n": 10000, "min_labels_before": 0,
                                                          "max_bytes": spage}, paging="all"), OWNER_LIBRARY_ONLY))
    if p.rich and cpage:
        p.add(lib(_op("claims-paged", "graphlet_claims", {"arm": a, "min_bp": 0, "route_consistent": False,
                                                          "labels": None, "max_bytes": cpage}, paging="all"),
                  OWNER_LIBRARY_ONLY))
    p.add(lib(_op("splits-n0", "graphlet_splits", {"arm": a, "n": 0, "min_labels_before": 0, "max_bytes": 2048}),
              OWNER_LIBRARY_ONLY))
    p.add(lib(_op("claims-unknown", "graphlet_claims", {"arm": a, "min_bp": 0, "route_consistent": True,
                                                        "labels": [{"name": "no such label"}], "max_bytes": 2048}),
              OWNER_LIBRARY_ONLY))
    p.add(lib(_op("continue-clash", "traverse_continue",
                  {"arm": a, "walk": ref("walks-length-%s" % a, "/rows/0/walk"), "overrides": {"arm": "right"},
                   "execute": False, "allow_unverified_index": False, "max_bytes": 16384, "reset_branches": False}),
              OWNER_NOT_SERVED))
    p.add(lib(_op("continue-oor", "traverse_continue",
                  {"arm": a, "walk": ref("walks-id-%s" % a, "/total"), "overrides": None, "execute": False,
                   "allow_unverified_index": False, "max_bytes": 16384, "reset_branches": False}), OWNER_NOT_SERVED))
    coords = bool(p.facts.get("coordinates"))
    p.add(lib(_op("walks-coords-%s" % a, "graphlet_walks", W(a, coordinates=True), paging="all"),
              REASON_COORDINATES))
    p.add(lib(_op("labels-coords-noname", "graphlet_labels", L(coordinates=True)), REASON_COORDINATES))
    if coords:
        p.add(lib(_op("labels-name-coords", "graphlet_labels",
                      L(arm=a, name=sel("top_label", labels="labels-%s-direct" % a), coordinates=True),
                      paging="all"), REASON_COORDINATES))
    p.add(lib(_op("walks-budget", "graphlet_walks", W(a, budget={"work_units": 1000})), REASON_BUDGET))
    if partner:
        p.add(lib(_op("exp-json-compare", "graphlet_export",
                      dict(E(format="json", what="compare"), path=None, other={"graphlet": partner})), REASON_EXPORT))


def plan_merges(p: Plan) -> None:
    """The displayed-parent cases: a label that entered a walk at a merge, filtered by route consistency."""
    for a in p.arms():
        if not any(o["id"] == "walk-merge-%s" % a for o in p.ops) or (p.facts.get("entered") or {}).get(a) is False:
            continue   # (entered False: the rehearsal saw every label entering the merge leave its first parent too)
        entered = sel("label_entered_at_merge", walk="walk-merge-%s" % a)
        for rc in (True, False):
            p.add(_op("walks-entered-%s-%s" % (a, "rc" if rc else "all"), "graphlet_walks",
                      W(a, label=entered, route_consistent=rc), paging="all"))
        p.add(_op("labels-entered-%s" % a, "graphlet_labels", L(arm=a, name=entered), paging="all"))
    for a in p.arms():
        top = sel("top_label", labels="labels-%s-direct" % a)
        for rc in (True, False):
            p.add(_op("walks-label-%s-%s" % (a, "rc" if rc else "all"), "graphlet_walks",
                      W(a, label=top, route_consistent=rc), paging="all"))


def plan_view_and_files(p: Plan) -> None:
    """graphlet_subtrie and the tools on its view; save / load; the list / free smoke (library-only)."""
    sid = p.sid
    a = p.arms()[0]
    top = sel("top_label", labels="labels-all-direct", index=0)
    p.add(lib(_op("sub-any", "graphlet_subtrie", {"labels": [top], "arm": None, "mode": "any"}), OWNER_NOT_SERVED))
    if p.labels() >= 2:
        p.add(lib(_op("sub-all", "graphlet_subtrie",
                      {"labels": [top, sel("top_label", labels="labels-all-direct", index=1)], "arm": a,
                       "mode": "all"}), OWNER_NOT_SERVED))
    view = ref("sub-any", "/handle")
    p.add(lib(_op("view-summary", "graphlet_summary", {"handle": view, "arm": None, "detail": "summary",
                                                       "max_bytes": MAXB}), OWNER_NOT_SERVED))
    p.add(lib(_op("view-walks", "graphlet_walks", dict(W(a), handle=view), paging="all"), OWNER_NOT_SERVED))
    p.add(lib(_op("view-walk", "graphlet_walk", {"handle": view, "arm": a, "walk": ref("view-walks", "/rows/0/walk"),
                                                 "max_bytes": MAXB}), OWNER_NOT_SERVED))
    p.add(lib(_op("view-labels", "graphlet_labels", dict(L(), handle=view), paging="all"), OWNER_NOT_SERVED))
    p.add(lib(_op("view-claims", "graphlet_claims", {"handle": view, "arm": a, "min_bp": 0, "route_consistent": True,
                                                     "labels": None, "max_bytes": MAXB}, paging="all"),
              OWNER_NOT_SERVED))
    p.add(lib(_op("view-splits", "graphlet_splits", {"handle": view, "arm": a, "n": 10000, "min_labels_before": 0,
                                                     "max_bytes": MAXB}, paging="all"), OWNER_NOT_SERVED))
    p.add(lib(_op("view-seq", "graphlet_sequence", dict(S(a, walk=ref("view-walks", "/rows/0/walk")), handle=view),
                  paging="all"), OWNER_NOT_SERVED))
    p.add(lib(_op("view-exp-fasta", "graphlet_export", dict(E(), path=None, other=None, handle=view)),
              OWNER_NOT_SERVED))
    p.add(lib(_op("view-exp-mgt", "graphlet_export", dict(E(format="mgt"), path=None, other=None, handle=view)),
              OWNER_NOT_SERVED))
    p.add(lib(_op("view-exp-gfa", "graphlet_export", dict(E(format="gfa"), path=None, other=None, handle=view)),
              OWNER_NOT_SERVED))
    p.add(lib(_op("view-compare", "graphlet_compare", dict(C(sid), a=view)), OWNER_NOT_SERVED))
    p.add(lib(_op("view-subtrie", "graphlet_subtrie", {"handle": view, "labels": [top], "arm": None, "mode": "any"}),
              OWNER_NOT_SERVED))
    p.add(lib(_op("sub-bad", "graphlet_subtrie", {"labels": [], "arm": None, "mode": "any"}), OWNER_NOT_SERVED))
    saved = "%s.saved.mgt" % sid
    p.add(lib(_op("save", "graphlet_save", {"path": saved}), OWNER_NOT_SERVED))
    p.add(lib(_op("load", "graphlet_load", {"path": saved}), OWNER_NOT_SERVED))
    loaded = ref("load", "/handle")
    p.add(lib(_op("load-summary", "graphlet_summary", {"handle": loaded, "arm": None, "detail": "summary",
                                                       "max_bytes": MAXB}), OWNER_NOT_SERVED))
    p.add(lib(_op("load-walks", "graphlet_walks", dict(W(a), handle=loaded), paging="all"), OWNER_NOT_SERVED))
    p.add(lib(_op("save-bad", "graphlet_save", {"path": "../saved.mgt"}), OWNER_NOT_SERVED))
    p.add(lib(_op("load-bad", "graphlet_load", {"path": "no-such-file.mgt"}), OWNER_NOT_SERVED))
    p.add(lib(_op("exp-nonmgt", "graphlet_export", dict(E(), path="%s.nonmgt.fa" % sid, other=None)),
              OWNER_NOT_SERVED))
    p.add(lib(_op("load-nonmgt", "graphlet_load", {"path": "%s.nonmgt.fa" % sid}), OWNER_NOT_SERVED))
    p.add(_op("smoke-list", "graphlet_list", {"max_bytes": MAXB}, ways=SMOKE_WAYS, not_exposed={"W3": OWNER_SMOKE}))
    p.add(_op("smoke-free", "graphlet_free", {"handle": loaded}, ways=SMOKE_WAYS, not_exposed={"W3": OWNER_SMOKE}))


def plan_compare(p: Plan, partner: str, full: bool, oid_prefix: str = "cmp") -> None:
    """The four modes against |partner| (a scenario id; the scenario's own id: a self-compare), with arm None and,
    when |full|, each arm; a bad mode."""
    arms = [None] + (p.arms() if full else [])
    for mode in ("claims", "walks", "labels", "prefix_subset"):
        for arm in arms:
            p.add(_op("%s-%s-%s" % (oid_prefix, mode, arm or "all"), "graphlet_compare", C(partner, arm=arm, mode=mode),
                      paging="all"))
    p.add(_op("%s-bad-mode" % oid_prefix, "graphlet_compare", C(partner, mode="all")))


REASON_FAILED_NO_DOC = "the service keeps no document for a failed seed"
FAILED_SEED_WAYS = ["W1", "W2", "W3", "W5"]


def plan_failed(p: Plan, other: str) -> None:
    """A failed seed: every tool the service exposes answers no_graphlet (W3, and W1 / W2 / W5: the runner's
    no_graphlet), one call each. W4 is not exposed (VERIFY-D VD2): the service keeps no document for a failed seed,
    so W4's runner would record no_document (a runner code: the op inconclusive by construction)."""
    for o in (_op("sum", "graphlet_summary", {"arm": None, "detail": "summary", "max_bytes": MAXB}),
              _op("walks", "graphlet_walks", W("left")),
              _op("walk", "graphlet_walk", {"arm": "left", "walk": 0, "max_bytes": MAXB}),
              _op("labels", "graphlet_labels", L()),
              _op("seq", "graphlet_sequence", S("left", walk=0)),
              _op("exp-mgt", "graphlet_export", E(format="mgt")),
              _op("cmp", "graphlet_compare", C(other))):
        o["ways"] = list(FAILED_SEED_WAYS)
        o["not_exposed"] = {"W4": REASON_FAILED_NO_DOC}
        p.add(o)


BIG_BODY_CASE = "graphlet.body_over_in_service_ceiling"   # async's case (service_by_design_g.json, 2026-10-20)
BIG_BODY_SUMMARY_CASE = BIG_BODY_CASE + ".summary"         # graphlet_summary answered from the stored row
MGT_EXPORT_CASE = "graphlet.export_mgt_is_the_stored_body"  # async's 4th case: W3's mgt receipt is the body's link
IN_SERVICE_MAX_BODY_BYTES = 8 * 1024 * 1024
LIBRARY_MAX_BYTES = 16 * 1024 * 1024                       # the cap case's upper bound (the library's own limit)
# the tools each case covers (async's BY-DESIGN-G-CASES.json graphlet.tools)
CAP_CASE_TOOLS = ("graphlet_summary", "graphlet_walks", "graphlet_walk", "graphlet_labels", "graphlet_compare")
BIG_BODY_TOOLS = ("graphlet_walks", "graphlet_walk", "graphlet_sequence", "graphlet_labels", "graphlet_compare",
                  "graphlet_export")


def by_design_case(o: dict, facts: dict) -> Optional[str]:
    """The by-design case an op W3 runs falls under (async's four area G cases), None when none applies:
    graphlet_export with format mgt (no arm, walks or what) the mgt case on any body (W3 answers the stored body's
    presigned link with its own receipt; the file is still compared by content); max_bytes above the service's
    262,144 (the cap's refusal comes first: the route refuses it before the handle, the token or the body is looked
    at); else, on a body over the in-service ceiling (facts over_8mib), graphlet_summary the .summary case (answered
    from the row) and the refusing tools the body case (graphlet_export only for fasta / gfa / json)."""
    ways = o.get("ways") or ALL_WAYS
    if "W3" not in ways or facts.get("failed"):
        return None   # a failed seed: every way answers no_graphlet, no case's sides can hold (VERIFY-D VD3)
    args = o["args"]
    if o["tool"] == "graphlet_export" and args.get("format") == "mgt":
        return MGT_EXPORT_CASE if all(args.get(k) is None for k in ("arm", "walks", "what", "export_id")) else None
    mb = o["args"].get("max_bytes")
    if isinstance(mb, bool) or not isinstance(mb, (int, type(None))):
        return None                    # a selector or a wrong type: no case's args match
    if mb is not None and MAXB < mb <= LIBRARY_MAX_BYTES and o["tool"] in CAP_CASE_TOOLS:
        return SERVICE_CAP_CASE
    if mb is not None and not 64 <= mb <= MAXB:
        return None
    if facts.get("over_8mib"):
        if o["tool"] == "graphlet_summary":
            return BIG_BODY_SUMMARY_CASE
        if o["tool"] in BIG_BODY_TOOLS and not (o["tool"] == "graphlet_export" and o["args"].get("format") == "mgt"):
            return BIG_BODY_CASE
    return None


def link_by_design(ops: List[dict], facts: dict) -> None:
    """Every op names the by-design case that applies to it (by_design_case), and no other."""
    for o in ops:
        if o.get("probe"):
            continue
        cid = by_design_case(o, facts)
        if cid:
            o["by_design"] = cid
        else:
            o.pop("by_design", None)


def plan_big(p: Plan) -> None:
    """A body over 8 MiB (DESIGN-v8 §4, D-R7: a six-call by-design plan): the service answers graphlet_summary from
    the stored record (in_service false, download) and every other tool but graphlet_export(mgt) body_too_large,
    while the library ways answer: each call a by-design pair (BIG_BODY_CASE), the export the body itself."""
    a = p.arms()[0]
    # the summary sends max_bytes 262,144: the row's stored summary next to the presigned link may not fit 2,048
    # (async's .summary case: else summary_trimmed)
    for o in (_op("sum", "graphlet_summary", {"arm": None, "detail": "summary", "max_bytes": MAXB}),
              _op("walks-id-%s" % a, "graphlet_walks", W(a), paging="all"),
              _op("walk-0-%s" % a, "graphlet_walk", {"arm": a, "walk": ref("walks-id-%s" % a, "/rows/0/walk"),
                                                     "max_bytes": MAXB}),
              _op("labels-all-direct", "graphlet_labels", L(), paging="all"),
              _op("seq-0-%s" % a, "graphlet_sequence", S(a, walk=ref("walks-id-%s" % a, "/rows/0/walk")),
                  paging="all")):
        p.add(o)                       # the cases: link_by_design (summary: .summary; the others: the body case)
    p.add(_op("exp-mgt", "graphlet_export", E(format="mgt")))


def plan_partner(p: Plan, a_id: str, full: bool, extras: dict) -> None:
    """A compare partner (DESIGN-v8 §4 pairs): its summary, the listings the compare's selectors and report read,
    and the compares against the corpus graphlet |a_id| (which runs first)."""
    p.add(_op("sum", "graphlet_summary", {"arm": None, "detail": "summary", "max_bytes": MAXB}))
    for a in p.arms():
        p.add(_op("walks-id-%s" % a, "graphlet_walks", W(a), paging="all"))
    p.add(_op("labels-all-direct", "graphlet_labels", L(), paging="all"))
    if not extras.get("no_compare"):   # the constrain / annotate partner: its compares run on the annotate graphlet
        plan_compare(p, a_id, full, "cmp")
    add_extras(p, extras)


def op_plan(sid: str, slot: str, facts: dict, rich: bool, pairs: Optional[List[Tuple[str, bool]]] = None,
            direct_only: bool = False, extras: Optional[dict] = None, partner_of: Optional[str] = None,
            reason: str = REASON_DIRECT_ONLY) -> List[dict]:
    """The op plan of one graphlet (DESIGN-v8 §4, INVENTORY §3.1): every distinct behaviour per tool and arm where it
    applies, per-walk sampling and paging passes only when |rich|; compares against |pairs| [(partner id, full)];
    |extras| the corpus-wide single ops this graphlet carries (cut kinds, probes, the by-design max_bytes call, the
    library-only view / files / smoke block). A direct-only graphlet runs W1, W2, W5 only; a compare partner
    (|partner_of|: the corpus graphlet it compares with) its own short plan."""
    extras = extras or {}
    p = Plan(sid, facts, rich)
    if partner_of:
        plan_partner(p, partner_of, bool(pairs and pairs[0][1]), extras)
    elif facts.get("failed"):
        plan_failed(p, extras.get("other") or sid)
    elif facts.get("over_8mib"):
        plan_big(p)
    else:
        plan_summary(p)
        plan_walks_and_walk(p)
        plan_labels(p)
        if facts.get("bases", True):
            plan_sequence(p)
        else:
            a = p.arms()[0]
            p.add(_op("seq-no-bases", "graphlet_sequence", S(a, walk=ref("walks-id-%s" % a, "/rows/0/walk"))))
        plan_merges(p)
        plan_export(p)
        plan_library_only(p, sid)
        for partner, full in pairs or []:
            plan_compare(p, partner, full, "cmp" if partner != sid else "self")
        add_extras(p, extras)
    link_by_design(p.ops, facts)
    if direct_only:
        kept = []
        for o in p.ops:
            if o.get("probe") or o.get("by_design") not in (None, MGT_EXPORT_CASE):
                continue   # W3's own behaviour: nothing to probe without an MCP job
            o.pop("by_design", None)   # the mgt export runs on the direct ways like any op; no W3, no case
            o["ways"] = [w for w in o.get("ways") or DIRECT_ONLY_WAYS if w in DIRECT_ONLY_WAYS]
            o["not_exposed"] = {"W3": reason, "W4": reason}
            kept.append(o)
        return kept
    return p.ops


def add_extras(p: Plan, extras: dict) -> None:
    """The corpus-wide single ops (one per kind in the whole corpus, DESIGN-v8 §5): deliberate cuts with cut_of, the
    W3-only probes, the by-design max_bytes 262,145 call, the library-only view / save / load / smoke block."""
    a = p.arms()[0]
    if extras.get("summary_cut"):
        cut = int(p.facts.get("summary_cut_bytes") or SUMMARY_CUT_BYTES)
        p.add(_op("sum-cut", "graphlet_summary", {"arm": None, "detail": "summary", "max_bytes": cut}, cut_of="sum"))
    paged = next((o for o in p.ops if o["id"] == "walks-paged-%s" % a), None)
    if extras.get("page_cut") and paged is not None:
        # the first page of a listing is the page cut (cut_of); its cursor reused with other arguments: bad_cursor
        page = dict(PAGED_WALKS, max_bytes=paged["args"]["max_bytes"])
        p.add(_op("walks-p1", "graphlet_walks", W(a, **page), paging={"pages": 1},
                  cut_of="walks-paged-%s" % a))
        p.add(_op("walks-reuse", "graphlet_walks", W(a, **dict(page, rank="length",
                                                               cursor={"cursor_of": "walks-p1"}))))
    if extras.get("ceiling_cut"):
        p.add(_op("walks-too-large", "graphlet_walks", W(a, max_bytes=64), cut_of="walks-id-%s" % a))
    if extras.get("row_cut"):
        p.add(_op("walks-full-cut", "graphlet_walks", W(a, rank="length", spell="full", max_bytes=ROW_PAGE),
                  paging="all", cut_of="walks-full-%s" % a))
    if extras.get("sequence_cut"):
        w = ref("walks-length-%s" % a, "/rows/0/walk")
        p.add(_op("seq-long", "graphlet_sequence", S(a, walk=w, with_seed=True), paging="all"))
        p.add(_op("seq-long-cut", "graphlet_sequence", S(a, walk=w, with_seed=True), cut_of="seq-long"))
    if extras.get("compare_cut"):
        partner = extras["compare_cut"]
        p.add(_op("cmp-cut", "graphlet_compare", C(partner, mode="claims", max_bytes=1024), paging="all",
                  cut_of="cmp-claims-all"))
    if extras.get("compare_paged"):
        cmp_page = page_of(p.facts, "cmp-paged", COMPARE_PAGE)
        if cmp_page:
            p.add(_op("cmp-paged", "graphlet_compare", C(extras["compare_paged"], mode="walks", max_bytes=cmp_page),
                      paging="all"))
    if extras.get("compare_paged_w3"):
        # a compare with rows that W3 pages (the rows pair, VERIFY-D-R4 R4-1): the service's compare paging path
        page = page_of(p.facts, "cmp-paged-w3", ROWS_PAIR_PAGE)
        if page:
            p.add(_op("cmp-paged-w3", "graphlet_compare", C(extras["compare_paged_w3"], mode="walks", max_bytes=page),
                      paging="all"))
    if extras.get("probes"):
        for probe in ("wrong_token", "bad_handle"):
            p.add(_op("probe-%s" % probe.replace("_", "-"), "graphlet_summary",
                      {"arm": None, "detail": "summary", "max_bytes": 2048}, ways=["W3"], probe=probe))
    if extras.get("service_cap"):
        p.add(_op("walks-cap", "graphlet_walks", W(a, max_bytes=MAXB + 1), paging="all", ways=["W1", "W2", "W3", "W4"],
                  by_design=SERVICE_CAP_CASE))
    if extras.get("self_compare") and not any(o["id"].startswith("self-") for o in p.ops):
        plan_compare(p, p.sid, False, "self")
    if extras.get("view_and_files"):
        plan_view_and_files(p)


def set_timeouts(ops: List[dict], timeouts: Optional[dict] = None) -> None:
    """Each op's timeout_s (FORMAT v2.3): the rehearsal's (rehearse_g: min(lane, max(floor, factor x median))) where
    it measured one, else the op's lane (DESIGN-v8 §3: 20 s interactive, 120 s heavy; graphlet_labels with a name runs
    on the heavy lane)."""
    for o in ops:
        t = (timeouts or {}).get(o["id"])
        if not isinstance(t, (int, float)) or isinstance(t, bool) or t <= 0:
            t = LANE_S["heavy"] if lane_of(o) == "heavy" else LANE_S["interactive"]
        o["timeout_s"] = float(t)


def lane_of(o: dict) -> str:
    if o["tool"] in HEAVY_TOOLS or (o["tool"] == "graphlet_labels" and o["args"].get("name") is not None):
        return "heavy"
    return "interactive"


# ===================================================================================================== scenarios
def default_corpus() -> dict:
    """The provisional corpus (no rehearsal yet): the first candidate of every slot (direct-only ones by default),
    conditional and v9 slots out; preconditions pending."""
    return {"version": 1, "provisional": True, "slots": {}, "candidates": {}, "pairs": {}}


def load_corpus(path: Optional[str]) -> dict:
    if not path:
        return default_corpus()
    with open(path) as f:
        c = json.load(f)
    if not isinstance(c, dict) or c.get("version") != 1:
        raise SystemExit("%s: not an area G corpus file (version 1)" % path)
    return c


def choose(cands: List[dict], corpus: dict, include_direct_only: bool = True) -> Dict[str, Optional[dict]]:
    """{slot: the candidate chosen, or None}: the corpus's choice, else (provisional) the first candidate of the slot
    (a direct-only one unless |include_direct_only| is false: the owner's decision of 2026-10-19 makes them v8's),
    conditional and v9 slots none."""
    by_slot: Dict[str, List[dict]] = {}
    for c in cands:
        by_slot.setdefault(c["slot"], []).append(c)
    out: Dict[str, Optional[dict]] = {}
    for slot, _stem, _title, conditional in SLOTS:
        pick = None
        if slot in V9_SLOTS:
            out[slot] = None
            continue
        stated = (corpus.get("slots") or {}).get(slot)
        if isinstance(stated, dict) and "candidate" in stated:
            pick = next((c for c in by_slot.get(slot, []) if c["cid"] == stated["candidate"]), None)
            if stated["candidate"] is not None and pick is None:
                raise SystemExit("corpus: slot %s names candidate %s, which the generator does not build" % (
                    slot, stated["candidate"]))
            if pick is not None and not expressible(pick) and not include_direct_only:
                pick = None   # a direct-only graphlet enters only with the owner's decision (--include-direct-only)
        elif not corpus.get("slots") or corpus.get("provisional"):
            if not conditional:
                pick = next((c for c in by_slot.get(slot, []) if include_direct_only or expressible(c)), None)
        out[slot] = pick
    return out


def facts_of(c: dict, corpus: dict) -> dict:
    """The facts the op plan reads: the rehearsal's measured ones, else the candidate's expected ones (v7's
    evidence)."""
    measured = ((corpus.get("candidates") or {}).get(c["cid"]) or {}).get("facts")
    f = dict(c["expected"])
    if isinstance(measured, dict):
        f.update(measured)
    if c["extra"].get("output", {}).get("sequences") is False:
        f["bases"] = False
    if c["extra"].get("output", {}).get("coordinates"):
        f["coordinates"] = True
    f.setdefault("radius", (c["overrides"].get("bounds") or {}).get("max_extension_bp"))
    sf = ((corpus.get("candidates") or {}).get(c["cid"]) or {}).get("service_facts") or {}
    if sf.get("body_over_8mib") or (f.get("body_bytes") or 0) > IN_SERVICE_MAX_BODY_BYTES:
        f["over_8mib"] = True
    if not f.get("arms"):
        d = c["overrides"].get("direction") or "both"
        f["arms"] = ["left", "right"] if d == "both" else [d]
    return f


def g_scenario(sid: str, title: str, c: dict, corpus: dict, library_commit: Optional[str], ops: List[dict],
               tags: List[str], timing_class: str = "fast") -> dict:
    """One area G scenario of candidate |c| (module docstring)."""
    measured = (corpus.get("candidates") or {}).get(c["cid"]) or {}
    strategy = candidate_strategy(c)
    echo = measured.get("strategy")
    if isinstance(echo, dict) and expressible(c):
        if echo != strategy:
            raise SystemExit("%s: the rehearsal's strategy echo differs from the merged strategy of %s: rebuild the "
                             "candidate from the echo (rehearsal report)" % (sid, c["cid"]))
    request = {"seeds": [copy.deepcopy(c["seed"])], "strategy": strategy}
    twin = twin_strategy(strategy)
    failed = bool(c["expected"].get("failed"))
    if failed:
        pre: dict = {"no_graphlet": True}
    else:
        pre = copy.deepcopy(measured.get("preconditions") or {})
    expect = [{"check": "status", "args": {"http": 200}}, {"check": "invariants", "args": {}},
              {"check": "traverse_seed_outcome", "args": {"walks": ["failed"] if failed else ["complete", "partial"]}}]
    set_timeouts(ops, measured.get("timeouts"))
    g: Dict[str, Any] = {
        "slot": c["slot"], "candidate": c["cid"], "twin": {"strategy": twin}, "preconditions": pre,
        "library_commit": library_commit, "ops": ops,
        "ways": {"default": list(ALL_WAYS if expressible(c) else DIRECT_ONLY_WAYS)},
        "not_exposed": dict(NOT_EXPOSED_TABLE),
    }
    if measured.get("service_facts"):
        g["facts"] = copy.deepcopy(measured["service_facts"])
    if measured.get("reference"):
        g["reference"] = copy.deepcopy(measured["reference"])
    if not expressible(c):
        g["not_exposed"]["*"] = direct_reason(c)
    truth = dict(c["truth"], candidate={"cid": c["cid"], "source": c["source"], "evidence": c["evidence"]})
    mcp = mcp_traverse(request)
    if c.get("direct_only_reason") and "tool" in mcp:
        mcp = {"not_implemented": c["direct_only_reason"]}   # an expressible strategy the corpus keeps direct-only
    tags = sorted(set(tags + ["graphlet", "area-g", "slot-" + c["slot"]] +
                      ([] if expressible(c) else ["direct-only"]) + (["conditional"] if c["conditional"] else []) +
                      ([] if measured.get("preconditions") or failed else ["rehearsal-pending"])))
    out = scenario(sid, AREA, title, "POST /traverse", request, truth, expect, timing_class=timing_class, tags=tags,
                   mcp=mcp, graphlet=g)
    if failed:
        # a failed seed walks nothing: no check that needs a walked seed (budget_respected, which scenario() adds for
        # a stated time budget, is inconclusive there by construction and blocks the gate: the shadow run)
        out["expect"] = [e for e in out["expect"] if e["check"] not in SEEDLESS_INCONCLUSIVE_CHECKS]
    return out


# the checks that need a walked seed (inconclusive on a failed seed: "no walked seed: nothing the budget governs")
SEEDLESS_INCONCLUSIVE_CHECKS = ("budget_respected",)


def gen_G(ctx: Any, built: List[dict], corpus: Optional[dict] = None, library_commit: Optional[str] = None,
          include_direct_only: bool = True) -> List[dict]:
    """The area G scenarios (module docstring): one per corpus slot, then the compare partners; ids G-<stem>-0001.
    |built| the scenarios built so far (v7's), whose T seeds the candidates reuse."""
    corpus = corpus if corpus is not None else default_corpus()
    cands = candidates(ctx, built)
    picks = choose(cands, corpus, include_direct_only)
    rich = set(corpus.get("rich") or DEFAULT_RICH)
    ids = dict((slot, "%s-0001" % stem) for slot, stem, _t, _c in SLOTS)
    present = dict((slot, c) for slot, c in picks.items() if c is not None)
    for slot, c in picks.items():
        if c is None:
            ctx.skip(AREA, ids[slot], "slot %s: %s" % (slot, V9_SLOTS.get(slot) or ((corpus.get("slots") or {}).get(
                slot) or {}).get("why") or ("conditional on the rehearsal" if dict((s[0], s[3]) for s in SLOTS)[slot]
                                            else "direct-only, left out (--no-direct-only)")))
    # the pairs: partner scenarios (built from the slot's candidate) and where the compare ops go
    pair_specs = pair_plan(present, ids, corpus)
    extras = corpus_extras(present, ids, pair_specs)
    if "failed_seed" in present:
        # a failed seed's compare names a graphlet that exists (its own resolves in no way: a selector error)
        graphlet = next((ids[s] for s in SLOT_NAMES if s in present and s != "failed_seed"), None)
        extras.setdefault(ids["failed_seed"], {})["other"] = graphlet
    out: List[dict] = []
    titles = dict((s[0], s[2]) for s in SLOTS)
    for slot, stem, _title, _cond in SLOTS:
        c = present.get(slot)
        if c is None:
            continue
        sid = ids[slot]
        facts = facts_of(c, corpus)
        pairs_here = [(ps["partner_of_a"], ps["full"]) for ps in pair_specs if ps["a"] == sid and ps["partner_of_a"]]
        ops = op_plan(sid, slot, facts, slot in rich, pairs_here, direct_only=not expressible(c),
                      extras=extras.get(sid), reason=direct_reason(c))
        out.append(g_scenario(sid, "area G %s: %s (%s)" % (slot.replace("_", " "), titles[slot], c["cid"]), c, corpus,
                              library_commit, ops, ["slot"], timing_class="slow" if slot in ("long_walk",
                                                                                              "over_8mib") else "fast"))
    for ps in pair_specs:
        if ps.get("partner_candidate") is None:
            continue
        pc = ps["partner_candidate"]
        facts = facts_of(pc, corpus)
        ops = op_plan(ps["partner_id"], pc["slot"], facts, False, [(ps["a"], ps["full"])],
                      direct_only=not expressible(pc), extras=extras.get(ps["partner_id"]), partner_of=ps["a"],
                      reason=direct_reason(pc))
        out.append(g_scenario(ps["partner_id"], "area G pair %s: %s (partner of %s)" % (ps["pair"], ps["title"],
                                                                                         ps["a"]),
                              pc, corpus, library_commit, ops, ["pair", "pair-" + ps["pair"]]))
    # partners first: a compare's b must have run (its MCP job) before the op that names it
    out = reorder_pairs(out, pair_specs)
    out += rows_pair(cands, corpus, library_commit)
    if "long_walk" in present:
        # the long walk's warm-ups right before it (the owner, 2026-10-19)
        i = next(j for j, s in enumerate(out) if s["id"] == ids["long_walk"])
        out[i:i] = long_walk_warm_ups(ctx, present["long_walk"], ids["long_walk"])
    return out


def _variant(c: dict, cid: str, overrides: Optional[dict] = None, seed: Optional[dict] = None,
             truth: Optional[dict] = None, source: str = "") -> dict:
    v = copy.deepcopy(c)
    v["cid"] = cid
    if overrides is not None:
        v["overrides"] = overrides
    if seed is not None:
        v["seed"] = seed
    if truth is not None:
        v["truth"] = truth
    v["source"] = source or v["source"]
    v["expected"] = dict((k, v2) for k, v2 in v["expected"].items() if k in ("arms", "labels", "failed"))
    return v


# A compare with rows on W3 (VERIFY-D-R4 R4-1): no MCP-expressible pair of the corpus differs within its common
# evidence (every compare of the shadow run counted 0 rows); of the rehearsed candidates on one seed, cut_lists.2 (a
# column walk of T-labels-0003's eco seed at 50 bp, cap 8) and over_8mib.1 (a trie of the same seed at 500 bp) do: on
# rehearsal 1's documents their walks compare counts 298 rows (claims 310, prefix_subset 300), all qualified, whole in
# 177,051 bytes; cmp-paged-w3 pages it (2 whole pages at 131,072, and at W3's budget less 512). Both run every way;
# b (the trie) first, without compares of its own.
ROWS_PAIR = {"a": "cut_lists.2", "b": "over_8mib.1", "id_a": "G-pair-rows-0002", "id_b": "G-pair-rows-0001"}
ROWS_PAIR_PAGE = 131072


def rows_pair(cands: List[dict], corpus: dict, library_commit: Optional[str]) -> List[dict]:
    """The rows pair's two scenarios (b, then a with its compares and cmp-paged-w3), [] when a candidate is missing
    or the corpus holds no rehearsal of it (a provisional corpus)."""
    by = dict((c["cid"], c) for c in cands)
    entries = corpus.get("candidates") or {}
    ca, cb = by.get(ROWS_PAIR["a"]), by.get(ROWS_PAIR["b"])
    if ca is None or cb is None or not (entries.get(ca["cid"]) or {}).get("preconditions") or \
            not (entries.get(cb["cid"]) or {}).get("preconditions"):
        return []
    ida, idb = ROWS_PAIR["id_a"], ROWS_PAIR["id_b"]
    title = "area G pair rows: a column walk at a small cap and a trie of one seed (compare rows, W3 paging)"
    fb = facts_of(cb, corpus)
    fb.pop("over_8mib", None)   # a pair partner, not the v9 slot
    ops_b = op_plan(idb, cb["slot"], fb, False, [(ida, False)], partner_of=ida, extras={"no_compare": True})
    ops_a = op_plan(ida, ca["slot"], facts_of(ca, corpus), False, [(idb, False)], partner_of=idb,
                    extras={"compare_paged_w3": idb})
    return [g_scenario(idb, title + " (b: %s)" % cb["cid"], cb, corpus, library_commit, ops_b, ["pair", "pair-rows"]),
            g_scenario(ida, title + " (a: %s)" % ca["cid"], ca, corpus, library_commit, ops_a, ["pair", "pair-rows"])]


def pair_plan(present: Dict[str, dict], ids: Dict[str, str], corpus: dict) -> List[dict]:
    """The compare pairs (PAIRS) whose slot is in the corpus: [{pair, title, a, partner_id, partner_candidate,
    partner_of_a, a_id_for_partner, full}]. The compare ops go on the partner (b = the corpus graphlet, which runs
    first) except the self-compare (on the graphlet itself)."""
    out = []
    for pair, stem, slot, title in PAIRS:
        c = present.get(slot)
        if c is None:
            continue
        a = ids[slot]
        if pair == "self":
            out.append({"pair": pair, "title": title, "a": a, "partner_id": None, "partner_candidate": None,
                        "partner_of_a": a, "a_id_for_partner": a, "full": False})
            continue
        if pair == "two_databases":
            continue   # conditional: the second database is named by the rehearsal (rehearse_g)
        if pair == "radius":
            ov = copy.deepcopy(c["overrides"])
            ov["bounds"]["max_extension_bp"] = max(31, ov["bounds"]["max_extension_bp"] // 2)
            pc = _variant(c, c["cid"] + "+radius", overrides=ov, source=c["source"] + " at half the radius")
        elif pair == "label_cap":
            ov = copy.deepcopy(c["overrides"])
            ov.setdefault("labels", {})["max_labels_per_node"] = 2
            pc = _variant(c, c["cid"] + "+cap2", overrides=ov,
                          source=c["source"] + " under labels.max_labels_per_node 2")
        elif pair == "constrain_annotate":
            pc = _variant(c, c["cid"] + "+constrain", source=c["source"] + " under labels.mode constrain (the walk "
                                                                          "preset's)")
            pc["extra"] = {}   # the walk preset as it is: constrain
            pc["direct_only_reason"] = REASON_ANNOTATE_PAIR   # a direct-only pair (VERIFY-D VD8, the coordinator)
            out.append({"pair": pair, "title": title, "a": a, "partner_id": "%s-0001" % stem,
                        "partner_candidate": pc, "partner_of_a": "%s-0001" % stem, "a_id_for_partner": a,
                        "full": False, "ops_on_a": True})
            continue
        elif pair == "two_seeds":
            t = _shift_seed(c["truth"], 10)
            if t is None:
                continue
            seed = dict(c["seed"], sequence=t["window"])
            pc = _variant(c, c["cid"] + "+shift10", seed=seed, truth=t, source=c["source"] + " shifted 10 bp right")
        else:
            continue
        out.append({"pair": pair, "title": title, "a": a, "partner_id": "%s-0001" % stem, "partner_candidate": pc,
                    "partner_of_a": None, "a_id_for_partner": a, "full": pair in ("radius", "two_seeds")})
    for ps in out:
        if ps["partner_id"] and not ps.get("ops_on_a"):
            ps["partner_of_a"] = None
    return out


def corpus_extras(present: Dict[str, dict], ids: Dict[str, str], pair_specs: List[dict]) -> Dict[str, dict]:
    """Where the corpus-wide single ops go (one per kind, DESIGN-v8 §5): {scenario id: extras}."""
    ex: Dict[str, dict] = {}

    def put(slot: str, **kw: Any) -> None:
        if slot in present:
            ex.setdefault(ids[slot], {}).update(kw)
    put("two_arms", summary_cut=True, ceiling_cut=True, page_cut=True, probes=True, service_cap=True)
    put("merges", view_and_files=True)
    if "long_walk" in present:
        put("long_walk", row_cut=True, sequence_cut=True)
    else:
        put("two_arms", row_cut=True)
    # the compare verdicts true and unknown (VERIFY-D VD8): self-compares of an uncut graphlet (merges.1: true), of a
    # depth-0 one and of one without bases (unknown)
    for slot in ("merges", "depth_zero", "no_bases"):
        put(slot, self_compare=True)
    paged = False
    for ps in pair_specs:
        if ps["pair"] == "constrain_annotate":
            # rows: constrain and annotate walk different paths (the claims compare pages at COMPARE_PAGE)
            ex.setdefault(ps["a"], {}).update(compare_paged=ps["partner_id"])
            ex.setdefault(ps["partner_id"], {}).update(no_compare=True)
            paged = True
    for ps in pair_specs:
        if ps["pair"] == "radius" and ps["partner_id"]:
            ex.setdefault(ps["partner_id"], {}).update(compare_cut=ps["a"])
            if not paged:
                ex[ps["partner_id"]]["compare_paged"] = ps["a"]
    return ex


def reorder_pairs(out: List[dict], pair_specs: List[dict]) -> List[dict]:
    """Slot scenarios in SLOTS order, each pair partner right after the graphlet it compares with (b runs first), or
    right before it when the compares run on the graphlet (ops_on_a: the constrain / annotate pair)."""
    partners = dict((ps["a"], []) for ps in pair_specs)
    before = dict((ps["a"], []) for ps in pair_specs)
    by_id = dict((s["id"], s) for s in out)
    for ps in pair_specs:
        if ps["partner_id"] and ps["partner_id"] in by_id:
            (before if ps.get("ops_on_a") else partners).setdefault(ps["a"], []).append(by_id[ps["partner_id"]])
    pids = set(ps["partner_id"] for ps in pair_specs if ps["partner_id"])
    res = []
    for s in out:
        if s["id"] in pids:
            continue
        res.extend(before.get(s["id"], []))
        res.append(s)
        res.extend(partners.get(s["id"], []))
    return res


# ===================================================================================================== expectations
# The library's answer to each error behaviour of the plans (INVENTORY §1.1, §2): the rehearsal's shadow run of W1
# (rehearse_g) and the tests hold every op to it, so an op that answers otherwise is named before the freeze. Keyed
# by the op id's stem (an id ending -left / -right / -all / a number keeps its stem).
EXPECTED_ERRORS = {
    "sum-bad-detail": "bad_argument", "sum-bad-arm": "bad_arm", "sum-arm-missing": "bad_arm",
    "walks-arm-none": "bad_arm", "walks-bad-rank": "bad_argument", "walks-n0": "bad_argument",
    "walks-unknown-label": "unknown_label", "walk-out": "bad_argument", "walk-neg": "bad_argument",
    "walks-reuse": "bad_cursor", "labels-unknown": "unknown_label", "seq-both": "bad_argument",
    "seq-neither": "bad_argument", "seq-end-lt-start": "bad_argument", "seq-start-past": "bad_argument",
    "seq-seed-seg": "bad_argument", "seq-seg-oor": "bad_argument", "seq-no-bases": "no_bases",
    "exp-bad-what": "bad_argument", "exp-path": "path_not_allowed", "splits-n0": "bad_argument",
    "claims-unknown": "unknown_label", "continue-clash": "bad_argument", "continue-oor": "bad_argument",
    "labels-coords-noname": "bad_argument", "walks-budget": "bad_argument", "view-compare": "view_unsupported",
    "view-subtrie": "view_unsupported", "view-exp-gfa": "view_unsupported", "sub-bad": "bad_argument",
    "save-bad": "path_not_allowed", "load-bad": "io_error", "load-nonmgt": "format_error",
    "cmp-bad-mode": "bad_argument", "self-bad-mode": "bad_argument", "walks-too-large": "result_too_large",
}


def expected_outcome(o: dict, facts: dict) -> str:
    """What the library answers to op |o| on a graphlet of |facts|: "ok" or an error code ("no_graphlet" for every
    op of a failed seed; walks with coordinates: no_coordinates unless the retrieval has coordinates)."""
    if facts.get("failed"):
        return "no_graphlet"
    oid = o["id"]
    if oid.startswith("walks-coords-"):
        return "ok" if facts.get("coordinates") else "no_coordinates"
    if facts.get("bases") is False and o["tool"] == "graphlet_export" and o["args"].get("format") in ("fasta", "gfa") \
            and oid not in EXPECTED_ERRORS:
        return "bad_argument"   # to_fasta / to_gfa need the bases (output.sequences: true)
    return EXPECTED_ERRORS.get(oid, "ok")


# ===================================================================================================== validation
def validate_g(scenarios: List[dict], cases: Optional[set] = None) -> Tuple[List[str], List[str]]:
    """(problems, pending) of the area G scenarios: problems refuse the build; pending lists what the freeze still
    needs (a by-design case async adds to the table, preconditions the rehearsal freezes, the library pin)."""
    problems: List[str] = []
    pending: List[str] = []
    ids = set(s["id"] for s in scenarios)
    for s in scenarios:
        g = s.get("graphlet")
        if not isinstance(g, dict):
            continue
        sid = s["id"]
        req = s["request"]
        st = req["strategy"]
        if st.get("output", {}).get("detail") != "graphlet":
            problems.append("%s: output.detail is not graphlet" % sid)
        if len(req.get("seeds") or []) != 1:
            problems.append("%s: not one seed" % sid)
        if (g.get("twin") or {}).get("strategy") != twin_strategy(st):
            problems.append("%s: the twin's strategy is not the request's with detail full and without "
                            "bounds.max_memory_mb" % sid)
        mcp = s.get("mcp") or {}
        if "tool" in mcp:
            a = mcp["arguments"]
            if mcp_merge(a.get("preset", "walk"), a.get("overrides")) != st:
                problems.append("%s: preset + overrides do not merge to request.strategy" % sid)
            if g["ways"]["default"] != ALL_WAYS:
                problems.append("%s: an MCP scenario runs the five ways" % sid)
        elif g["ways"]["default"] != DIRECT_ONLY_WAYS:
            problems.append("%s: a direct-only scenario runs W1, W2, W5" % sid)
        tb = (st.get("bounds") or {}).get("time_budget_ms")
        if not isinstance(tb, (int, float)):
            problems.append("%s: no time budget in the strategy" % sid)
        if not g.get("library_commit"):
            pending.append("%s: graphlet.library_commit (the pin: build_manifest.py --library-commit)" % sid)
        if not g.get("preconditions"):
            pending.append("%s: preconditions (the rehearsal freezes them)" % sid)
        seen = set()
        # over the in-service ceiling: the rehearsal's service facts, else the slot's expectation (facts_of)
        big = (g.get("facts") or {}).get("body_over_8mib")
        gf = dict(g, _facts={"over_8mib": bool(big if big is not None else g.get("slot") == "over_8mib"),
                             "failed": bool((g.get("preconditions") or {}).get("no_graphlet"))})
        for o in g.get("ops") or []:
            oid = o.get("id")
            if oid in seen:
                problems.append("%s: duplicate op %s" % (sid, oid))
            seen.add(oid)
            problems += ["%s/%s: %s" % (sid, oid, m) for m in op_problems(o, gf, ids, seen)]
            bd = o.get("by_design")
            if bd and cases is not None and bd not in cases:
                pending.append("%s/%s: by-design case %s not in the by-design tables yet (async adds it)" % (
                    sid, oid, bd))
    return problems, sorted(set(pending))


# An op's raw pages are raw/<scenario>.<op>.<page>.json.gz in its way's directory, and W3's way directory is the MCP
# arm's (arms/mcp), whose execution record is raw/<scenario>.mcp.<repeat>.json.gz (suite_common.raw_name): an op
# named after an arm would overwrite it (async's runner exits 4 on such a plan, 2026-10-20). Op ids keep to the
# characters a raw name keeps (suite_common.raw_name / run_graphlet_local.safe), so no two ops share a page name.
OP_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


def op_id_problems(oid: Any) -> List[str]:
    if not isinstance(oid, str) or not OP_ID.match(oid):
        return ["op id %r: letters, digits, '.', '_' and '-' only (its raw page names)" % (oid,)]
    if oid in sc.ARMS:
        return ["op id %r is an arm's name: its raw pages would collide with the execution's raw record "
                "(raw/<scenario>.%s.<repeat>.json.gz)" % (oid, oid)]
    return []


def op_problems(o: dict, g: dict, ids: set, before: set) -> List[str]:
    """What is wrong with one op: its tool and arguments for the ways it runs (the service's arguments all named on
    an op W3 runs, no library-only argument there; a library-only tool or argument without W3 and with the owner's
    reason), its paging, selectors naming earlier ops, a cut_of naming an earlier op, timeout_s."""
    out = []
    tool, args = o.get("tool"), o.get("args")
    ways = o.get("ways") or (g.get("ways") or {}).get("default") or ALL_WAYS
    out += op_id_problems(o.get("id"))
    if not isinstance(args, dict):
        return out + ["args is not an object"]
    if o.get("paging") not in ("none", "all") and not (isinstance(o.get("paging"), dict)
                                                      and set(o["paging"]) == {"pages"}):
        out.append("paging %r" % (o.get("paging"),))
    if "W3" in ways and not o.get("probe"):
        if tool not in SERVICE_ARGS:
            out.append("%s is not a service tool, yet W3 runs it" % tool)
        else:
            want = set(SERVICE_ARGS[tool])
            have = set(args) - {"cursor"}
            if have != want:
                out.append("the service's arguments %s, the op names %s" % (sorted(want), sorted(have)))
    if "W3" not in ways:
        reason = (o.get("not_exposed") or {}).get("W3") if isinstance(o.get("not_exposed"), dict) else None
        if not reason:
            out.append("W3 left out without the owner's reason")
    if tool in LIBRARY_TOOL_ARGS:
        if "W3" in ways:
            out.append("library-only tool %s on W3" % tool)
        want = set(LIBRARY_TOOL_ARGS[tool]) | ({"handle"} if tool != "graphlet_list" else set())
        if not set(args) <= want or not set(LIBRARY_TOOL_ARGS[tool]) <= set(args) | {"handle"}:
            out.append("the library's arguments %s, the op names %s" % (sorted(LIBRARY_TOOL_ARGS[tool]),
                                                                       sorted(args)))
    for k in LIBRARY_ONLY_ARGS.get(tool, ()):
        if k in args and "W3" in ways:
            out.append("library-only argument %s on W3" % k)
    for name in _selector_ops(args):
        if name not in before or name == o.get("id"):
            out.append("a selector reads op %r, which does not run before it" % name)
    if o.get("cut_of") and o["cut_of"] not in before:
        out.append("cut_of %r is not an earlier op" % o["cut_of"])
    if o.get("probe") and o.get("ways") != ["W3"]:
        out.append("a probe runs on W3 only")
    mb = args.get("max_bytes")
    if o.get("probe") and isinstance(mb, int) and mb > MAXB:
        # the service refuses max_bytes above 262,144 before it checks the token or the handle
        out.append("a probe sends max_bytes %d above %d: the cap's refusal would answer, not the probe's" % (mb, MAXB))
    if not o.get("probe"):
        want = by_design_case(o, g.get("_facts") or {})
        if (o.get("by_design") or None) != want:
            out.append("by_design %r, the case that applies %r" % (o.get("by_design"), want))
    t = o.get("timeout_s")
    if not isinstance(t, (int, float)) or isinstance(t, bool) or t <= 0:
        out.append("timeout_s %r" % (t,))
    return out


def _selector_ops(v: Any) -> List[str]:
    """The op ids the selectors in |v| read (ref, cursor_of, the named selectors' op arguments)."""
    out: List[str] = []
    if isinstance(v, dict):
        if set(v) == {"ref", "path"} and isinstance(v["ref"], str):
            out.append(v["ref"])
            return out
        if set(v) == {"cursor_of"}:
            out.append(v["cursor_of"])
            return out
        if "select" in v and set(v) <= {"select", "args"}:
            for k, x in (v.get("args") or {}).items():
                if k in ("walks", "walk", "labels", "of", "length_of") and isinstance(x, (str, list)):
                    out += [x] if isinstance(x, str) else [y for y in x if isinstance(y, str)]
                else:
                    out += _selector_ops(x)
            return out
        for x in v.values():
            out += _selector_ops(x)
    elif isinstance(v, list):
        for x in v:
            out += _selector_ops(x)
    return out


def summary(scenarios: List[dict]) -> dict:
    """Counts of area G: scenarios, ops, calls per way (an op = one call per way it runs; paging aside)."""
    per_way: Dict[str, int] = {}
    n_ops = 0
    for s in scenarios:
        g = s.get("graphlet")
        if not isinstance(g, dict):
            continue
        for o in g["ops"]:
            n_ops += 1
            for w in o.get("ways") or g["ways"]["default"]:
                per_way[w] = per_way.get(w, 0) + 1
    return {"scenarios": sum(1 for s in scenarios if isinstance(s.get("graphlet"), dict)), "ops": n_ops,
            "calls_per_way": dict(sorted(per_way.items()))}
