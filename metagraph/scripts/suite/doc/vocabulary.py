"""The one vocabulary of the results document (the readability review of the full v4 report, READ-022, READ-003,
READ-007, READ-008): the three layers an outcome has, each with its words, defined once and rendered from here by
the reading guide's legend (I.3), appendix F's tables and the glossary's entries (appendix D, through the `vocab`
placeholders), so that the three cannot drift. The definitions follow checks.evaluate (the check items),
suite_common.outcome_from_checks, the runners and reevaluate._outcome_of (the execution outcome), cells.arm_state
(the displayed standing) and compare.py (the classes); where the layers differ, the difference is stated, not
smoothed over: a check that raises inside the evaluation is a failing item (checks.evaluate catches it), so the
execution reads `fail` on every access path; only where the evaluation itself raised is the execution recorded
`error` on a direct access path and `fail` via MCP (run_direct.py, reevaluate.py).

Every text is plain prose (escaped by the renderer); `code` marks a code span for latex.inline.
"""

ARMS = ("direct-local", "direct-remote", "mcp")

# the three access paths and the one word the records use for them (READ-008)
ACCESS_PATH = ("The suite's records, the tables' row labels and the stored keys call the three access paths to the "
               "host *arms*: direct-local, direct-remote and mcp. In this document's prose they are access paths; "
               "a *traversal arm* is a direction of a walk (left or right of the seed) and is always written with "
               "that qualifier.")

# layer 1: the outcome of one check item (checks.evaluate; suite_common.CHECK_OUTCOMES)
CHECK_OUTCOMES = (
    ("pass", "the check tested its claim on this answer and the claim held"),
    ("fail", "the claim did not hold, or the check could not run as a check: an unknown check name or a check "
             "that raised is recorded as a failing item stating why, never a silent pass"),
    ("inconclusive", "the check could not test its claim on this answer (an allowed 503 deadline, a list a stated "
                     "hiding shortened, a walk a stated limit stopped first, nothing to compare)"),
    ("finding", "a stated bound was broken on the clock the check reads: a 200 past its budget (the host's own "
                "elapsed time), a 503 deadline more than max(1 s, the budget) late (the client's clock: the suite's "
                "tolerance), a walk stopped early or run past its pieces; whatever the HTTP status"),
    ("not_exposed", "the MCP answer lacks the field the check reads: recorded, neither passed nor failed, and not "
                    "compared (on a direct access path a missing field is a failure, the field being part of the "
                    "contract)"),
)

# layer 2: the outcome of one execution (suite_common.OUTCOMES), with how it is decided and what it means
EXECUTION_OUTCOMES = (
    {"key": "pass", "legend": "every check that ran passed and at least one tested the claim (a property check, or the generic checks of a test without one)",
     "glyph": "\\gpass", "word": "pass", "decided": "by the checks",
     "meaning": "every check that ran passed and at least one of them tested the claim: a property check, or, for a "
                "scenario without one, its generic checks (a refusal: its status, code and text are the claim; a "
                "warm-up: what every answer must hold); a "
                "check that did not apply to the answer, or that a property-level mark of the manifest set aside, "
                "is not counted",
     "consequence": "the property was shown on this execution"},
    {"key": "fail", "legend": 'a check failed: the claim did not hold, a required field was missing on a direct access path, or a check could not run as a check',
     "glyph": "\\gfail", "word": "fail", "decided": "by the checks",
     "meaning": "a check failed (its claim did not hold, a required field was missing on a direct access path, a "
                "stated stop or hiding could not have happened, or a check could not run as a check)",
     "consequence": "a functional failure: the verdict is FAIL and the gate fails"},
    {"key": "error", "legend": 'no HTTP answer (a transport failure, a client time-out, an undecodable body); the evaluation itself (checks.evaluate) raised on a direct access path; a recorded outcome outside the vocabulary',
     "glyph": "\\gerr", "word": "error", "decided": "by the runner, without the checks",
     "meaning": "no HTTP answer to check: a transport failure, a client time-out or an undecodable body; where the "
                "evaluation itself (checks.evaluate) raised on a direct access path the execution is recorded as an "
                "error too (via MCP the same line reads fail), while a single check that raises inside the evaluation "
                "is a failing item, so that execution reads fail on every access path; an outcome a runner recorded "
                "outside this vocabulary is read as an error with the recorded word kept",
     "consequence": "counted with the failures"},
    {"key": "finding", "legend": "a stated bound broken on the clock the check reads (a 200 past its budget by the host's own elapsed time; an allowed 503 deadline more than max(1 s, the budget) late on the client's clock, the suite's tolerance; a walk stopped early or run past its pieces): apart from fail, whatever the HTTP status",
     "glyph": "\\gfind", "word": "finding", "decided": "by the checks",
     "meaning": "a stated bound was broken on the clock the check reads: a 200 past its budget (the host's own "
                "elapsed time), an allowed 503 deadline more than max(1 s, the budget) late (the client's clock: the "
                "suite's tolerance for a deadline answer, not a proof of where the time went), a walk stopped before "
                "its budget or run past it beyond its pieces. The HTTP status may be a 503 that carried no results: "
                "finding is the suite's category for a resource-budget breach, apart from the functional-check "
                "failures",
     "consequence": "the verdict reads FINDINGS; on a gated access path it fails the gate"},
    {"key": "inconclusive", "legend": 'the answer could not show the property (an allowed 503 deadline without results, a shortened list, a walk a stated limit stopped, a stop that did not happen, nothing to compare, no field exposed via MCP, a host state): never a pass, never a failure',
     "glyph": "\\ginc", "word": "inconclusive", "decided": "by the checks",
     "meaning": "the answer could not show the property: an allowed 503 deadline (no results; the SPEC allows it "
                "under load where the request names a budget), a list a stated hiding shortened, a walk a stated "
                "limit stopped first, a stop that was the property and did not happen, nothing to compare, no "
                "property field exposed via MCP, or a host state (the service ended the call by its own limit). "
                "Also when no check could test anything",
     "consequence": "never a pass and never a failure; on a gated access path it blocks the gate"},
    {"key": "not_implemented_mcp", "legend": "not exercised via MCP: the tools cannot express the request (the manifest's mark or a table case): a known limit, not a failure and not coverage",
     "glyph": "\\gni", "word": "not implemented via MCP",
     "decided": "by the runner, without the checks",
     "meaning": "the request cannot be asked through the MCP tools: the manifest marks it, or a by-design case of "
                "that agreement names the missing argument; the service's own refusal, when one was recorded, is "
                "kept with the line",
     "consequence": "a known limit of the tools, not a failure and not coverage: the cell reads not exercised, "
                    "never ok"},
    {"key": "service_by_design", "legend": "the service answered or refused in the host's place as its table documents; u: the pair is unverified; nothing of the property was tested through the tools there",
     "glyph": "\\gbd", "word": "service-side by design",
     "decided": "by the runner, replayed from the by-design table",
     "meaning": "the service answered or refused in the host's place as its table documents (a refusal at submit, a "
                "relayed host refusal under the service's code, an answer where the host refuses, a clamped value); "
                "the pair is verified when the direct answer matched the case's host side too",
     "consequence": "a documented difference, not a failure; nothing of the property was tested through the tools "
                    "on that execution"},
)

# layer 3: what a cell shows (cells.arm_state)
STANDING = ("The glyph of a cell is the test's standing on that access path: the worst execution outcome over its "
            "repeats (repeats that differ carry the mark m), except that on the MCP path the agreement class takes "
            "the outcome's place where the service decided the call (not exercised via MCP, service-side by design). "
            "A standing count counts tests. I.4's outcome table counts executions; a chapter's outcomes row counts "
            "tests too, each test's outcome over its repeats as the comparison records it (the common outcome, else "
            "the worst of error, fail, finding and inconclusive; mixed where the repeats differ otherwise).")

DIRECT_CLASSES = (
    ("identical", "equal answers (the same digest)"),
    ("same_checks", "different answers, the same check outcomes"),
    ("timing", "the answers differ by timing alone (one time-limited)"),
    ("not_comparable", "no check both access paths evaluated"),
    ("disagree", "a check differs: fails the run"),
)

MCP_CLASSES = (
    ("agree", "the checks both access paths evaluated agree (the depth says how many property checks)"),
    ("not_implemented_mcp", "by the manifest or the table (not exercised)"),
    ("service_by_design", "a table case matched (verified when both sides match the case)"),
    ("not_comparable", "no check both access paths evaluated"),
    ("disagree", "a check differs: fails the run"),
    ("missing", "no MCP execution recorded"),
)

# the three classifications a reader meets on one answer (READ-003)
CLASSIFICATIONS = ("Three classifications are distinct and are never read for one another: the HTTP status of the "
                   "answer (a 200 with results, a 4xx refusal, a 503 deadline that carried no results), the suite's "
                   "outcome of the execution (what the checks could show of the claim: pass, fail, error, finding, "
                   "inconclusive, or a mark the runner recorded), and the release decision (the deploy gate's "
                   "judgement over a selected population).")


def execution_outcome(key):
    for e in EXECUTION_OUTCOMES:
        if e["key"] == key:
            return e
    return None


# A run whose manifest has an answer-side scenario-level mark (blocks.scenario_mark_parts: the request is expressible,
# the MCP answer lacks every field the test's property checks read; v6's P-scope-0004 / -0006) widens the definition
# of not implemented via MCP by this clause; a run without one keeps the definition as it was (the delta verifier's
# M2, 2026-10-18).
ANSWER_SIDE_LEGEND = ("not exercised via MCP: the tools cannot express the request (the manifest's mark or a table "
                      "case), or, by the manifest's scenario-level mark, the MCP answer lacks every field the test's "
                      "property checks read: a known limit, not a failure and not coverage")
ANSWER_SIDE_MEANING = ("the request cannot be asked through the MCP tools: the manifest marks it, or a by-design case of "
                       "that agreement names the missing argument; or, by the manifest's scenario-level mark, the "
                       "request is expressible but the MCP answer lacks every field the test's property checks read; "
                       "the service's own refusal, when one was recorded, is kept with the line")


def legend_rows(answer_side=False):
    """[(glyph macro, word, meaning)] for the reading guide's glyph legend: the execution outcomes in the tables'
    glyph order and the one glyph for no execution; `answer_side`: the run has an answer-side scenario-level mark."""
    rows = [(e["glyph"], e["word"], ANSWER_SIDE_LEGEND if (answer_side and e["key"] == "not_implemented_mcp")
             else e["legend"]) for e in EXECUTION_OUTCOMES]
    rows.append(("\\gnone", "no execution", "the plan had none (a long scenario on direct-remote); a planned "
                                              "execution not recorded is footnoted"))
    return rows


def outcome_table_rows(answer_side=False):
    """[(key, decided, meaning, consequence)] for appendix F's execution-outcome table (`answer_side` as for
    legend_rows)."""
    return [(e["key"], e["decided"], ANSWER_SIDE_MEANING if (answer_side and e["key"] == "not_implemented_mcp")
             else e["meaning"], e["consequence"]) for e in EXECUTION_OUTCOMES]


def _sentence(text):
    """A clause as a sentence: its first letter capitalised (the glossary joins a meaning and a consequence)."""
    return text[:1].upper() + text[1:] if text else text


def glossary_texts():
    """The glossary's shared definitions (appendix D), plain text, filled through the `vocab` placeholders."""
    outcomes = "; ".join("%s (%s)" % (e["word"], e["decided"]) for e in EXECUTION_OUTCOMES)
    return {
        "check_outcome": "What one check item answers: " + "; ".join(
            "%s: %s" % (k, m.split(" (", 1)[0].split(":", 1)[0]) for k, m in CHECK_OUTCOMES)
        + ". A finding or an inconclusive item never counts as a pass.",
        "execution_outcome": "From the check items by precedence fail, finding, inconclusive, pass; inconclusive also "
                             "when no check could test anything. Decided without the checks: error, not implemented "
                             "via MCP, service-side by design. The seven outcomes: " + outcomes + ".",
        "error": execution_outcome("error")["meaning"] + ". " + _sentence(execution_outcome("error")["consequence"]) + ".",
        "finding": execution_outcome("finding")["meaning"] + ". " + _sentence(execution_outcome("finding")["consequence"]) + ".",
        "inconclusive": execution_outcome("inconclusive")["meaning"] + ". "
                        + _sentence(execution_outcome("inconclusive")["consequence"]) + ".",
        "pass": execution_outcome("pass")["meaning"] + ".",
        "standing": STANDING,
        "classifications": CLASSIFICATIONS,
        "access_path": ACCESS_PATH,
    }
