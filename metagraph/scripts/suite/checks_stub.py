#!/usr/bin/env python3
"""Stand-in for lane S1's checks.py with the same API (DESIGN §8.2), evaluating the `status` check only.

run_direct.py imports checks.py when it exists and falls back to this module, so the runner, the comparator
and the report can be built and smoke-tested before the real checks land. `status` is evaluated as S1's
check_status does (http, a number or a list; then the error's code, text, text_prefix, text_contains where
stated); every other expected check is reported `not_exposed` with the detail "not evaluated: checks stub".
run.json names the module that ran, and run_direct refuses the stub for a selection that expects any other
check unless --allow-stub, so a stub run is never mistaken for a verified one.
"""

STUB = True


def _status(answer, args):
    got = answer.get("http")
    want = args.get("http", 200)
    wants = want if isinstance(want, list) else [want]
    err = answer.get("error") or {}
    text = err.get("text") or ""
    if got not in wants:
        return {"check": "status", "outcome": "fail",
                "detail": "http %s, expected %s (error %s: %s)" % (got, want, err.get("code"), text[:300])}
    problems = []
    if "code" in args and err.get("code") != args["code"]:
        problems.append("code %r, expected %r" % (err.get("code"), args["code"]))
    if "text" in args and text != args["text"]:
        problems.append("text %r, expected %r" % (text[:400], args["text"]))
    if "text_prefix" in args and not text.startswith(args["text_prefix"]):
        problems.append("text %r does not start with %r" % (text[:400], args["text_prefix"]))
    for frag in args.get("text_contains", []) or []:
        if frag not in text:
            problems.append("text %r lacks %r" % (text[:400], frag))
    if problems:
        return {"check": "status", "outcome": "fail", "detail": "; ".join(problems)}
    return {"check": "status", "outcome": "pass",
            "detail": "http %s%s" % (got, (" code %s" % err.get("code")) if err.get("code") else "")}


def evaluate(scenario, answer, arm, previous=None, client=None):
    """[{check, outcome: pass|fail|not_exposed, detail}] for the scenario's `expect` list.

    `previous` (the answer digests of the scenario's earlier repeats on this arm) is accepted as lane S1's
    checks.py takes it, and not used here.

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
