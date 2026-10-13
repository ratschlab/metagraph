#!/usr/bin/env python3
"""Stand-in for lane S1's checks.py with the same API (DESIGN §8.2), checking the HTTP status only.

run_direct.py imports checks.py when it exists and falls back to this module, so the runner, the comparator
and the report can be built and smoke-tested before the real checks land. Every expected check other than
`status` is reported `not_exposed` with the detail "not evaluated: checks stub", and run.json names the
module that ran, so a stub run is never mistaken for a verified one.
"""

STUB = True


def _status(answer, args):
    want = args.get("http", 200)
    got = answer.get("http")
    ok = got == want
    detail = "http %s, expected %s" % (got, want)
    if args.get("code") is not None:
        detail += " (error code not compared by the stub)"
    return {"check": "status", "outcome": "pass" if ok else "fail", "detail": detail}


def evaluate(scenario, answer, arm):
    """[{check, outcome: pass|fail|not_exposed, detail}] for the scenario's `expect` list.

    Without a `status` expectation the answer's status must be 2xx (an implicit check, named `status`).
    """
    results = []
    saw_status = False
    for exp in scenario.get("expect") or []:
        name = exp.get("check")
        args = exp.get("args") or {}
        if name == "status":
            saw_status = True
            results.append(_status(answer, args))
        else:
            results.append({"check": name, "outcome": "not_exposed",
                            "detail": "not evaluated: checks stub (HTTP status only)"})
    if not saw_status:
        got = answer.get("http")
        ok = isinstance(got, int) and 200 <= got < 300
        results.insert(0, {"check": "status", "outcome": "pass" if ok else "fail",
                           "detail": "http %s, expected 2xx (implicit: no status expectation)" % got})
    return results
