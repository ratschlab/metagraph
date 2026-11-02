#!/usr/bin/env python3
"""Check mechanically detectable writing rules with a ratchet.

The scanner reads the files that `scripts/text-rules.json` names, pulls out the text a
user reads (prose in Markdown and reStructuredText, strings in code that reach
a user), and counts rule hits per file. `text-rules-baseline.json` holds the
counts that existing text already has (`scripts/text-rules-baseline.json`). The check fails when a file has more
hits than its baseline, or fewer (lower the baseline in the same change).

Usage:
  check_text_rules.py                    # check against the baseline
  check_text_rules.py --report           # list every hit, then the counts
  check_text_rules.py --update-baseline  # record reductions only
  check_text_rules.py --report FILE...   # list hits in some files only

Hits that only count, never fail: R2 (sentences over the word limit) and
R6 (passive voice). Mark a deliberate exception with `text-rules-allow: <rule>`
on the same line or the line before (in a comment).
"""
import argparse
import ast
import fnmatch
import json
import pathlib
import re
import subprocess
import sys

RATCHET = ("R1", "R4", "R5", "R7", "R8", "R9")
REPORT_ONLY = ("R2", "R6")
WORD = re.compile(r"[A-Za-z][A-Za-z'-]*")
PROSE = re.compile(r"[A-Za-z]{2,},? [A-Za-z]{2,},? [A-Za-z]{2,}")  # 3 plain words in a row
ALLOW = re.compile(r"text-rules-allow:\s*([R0-9, ]+)")


# ---------- extractors: each yields (line, text) ----------

def _blank(m):
    return "\n" * m.group(0).count("\n")


def markdown_units(src, rst=False):
    """Paragraphs and list items, without code, URLs and markup."""
    src = re.sub(r"<!--.*?-->", _blank, src, flags=re.S)
    lines = src.split("\n")
    out, buf, start = [], [], 0
    fence = None
    literal_indent = None
    directive_indent = None

    def flush():
        nonlocal buf
        if buf:
            out.append((start, " ".join(buf)))
        buf = []

    for i, raw in enumerate(lines, 1):
        line = raw.rstrip()
        stripped = line.lstrip()
        indent = len(line) - len(stripped)
        if rst and directive_indent is not None:
            if stripped and indent <= directive_indent:
                directive_indent = None
            elif re.match(r":[\w-]+:", stripped):
                continue
        if fence:
            if stripped.startswith(fence):
                fence = None
            continue
        m = re.match(r"(```+|~~~+)", stripped)
        if m:
            flush()
            fence = m.group(1)
            continue
        if rst and literal_indent is not None:
            if not stripped or indent > literal_indent:
                continue
            literal_indent = None
        if rst and (re.match(r"\.\. (code-block|code|sourcecode|highlight|literalinclude|image|figure|toctree|raw|automodule|autoclass|autofunction)\b", stripped)
                    or (line.endswith("::") and not stripped.startswith(".. "))):
            text = line[:-2] if line.endswith("::") else ""
            if text.strip() and not stripped.startswith(".."):
                buf.append(text.strip())
                if not start:
                    start = i
            flush()
            literal_indent = indent
            continue
        if not stripped or re.match(r"^[=\-~^*#+|:` ]{3,}$", stripped) or stripped.startswith("|"):
            flush()
            continue
        if rst and stripped.startswith(".. "):
            flush()
            directive = re.match(r"\.\. (?:note|warning|attention|caution|danger|error|hint|important|tip)::\s*(.*)", stripped)
            if directive:
                directive_indent = indent
                if directive.group(1):
                    out.append((i, directive.group(1)))
            continue
        if not rst and indent >= 4 and not buf:
            continue  # indented code block
        if re.match(r"([-*+]|\d+[.)]|#+)\s", stripped) or stripped.startswith("<"):
            flush()
        if not buf:
            start = i
        buf.append(stripped)
    flush()
    cleaned = []
    for ln, text in out:
        text = re.sub(r"``.*?``|`[^`]*`", " CODE ", text)
        text = re.sub(r"\]\([^)]*\)|<https?://[^>]*>|https?://\S+", "]", text)
        text = re.sub(r"<[^>]+>", " ", text)
        text = re.sub(r"&[a-z]+;|&#\d+;", " ", text)
        text = re.sub(r"^([-*+]|\d+[.)]|#+)\s+", "", text)
        cleaned.append((ln, text))
    return cleaned


def python_units(src, cfg):
    """Strings passed to the calls and keywords named in the config."""
    try:
        tree = ast.parse(src)
    except SyntaxError:
        return []
    calls = set(cfg.get("python_calls", []))
    keywords = set(cfg.get("python_keywords", []))
    dict_keys = set(cfg.get("python_dict_keys", []))
    decorators = set(cfg.get("python_docstring_decorators", []))
    assigns = set(cfg.get("python_assign_names", []))
    out = []

    def text_of(node):
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            return node.value
        if isinstance(node, ast.JoinedStr):
            return "".join(v.value if isinstance(v, ast.Constant) else "X" for v in node.values)
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
            a, b = text_of(node.left), text_of(node.right)
            return (a or "X") + (b or "X") if (a or b) else None
        return None

    def name_of(func):
        if isinstance(func, ast.Name):
            return func.id
        if isinstance(func, ast.Attribute):
            return func.attr
        return None

    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            name = name_of(node.func)
            if name in calls:
                for a in node.args[:1]:
                    t = text_of(a)
                    if t:
                        out.append((a.lineno, t))
            for kw in node.keywords:
                if kw.arg in keywords:
                    t = text_of(kw.value)
                    if t:
                        out.append((kw.value.lineno, t))
        elif isinstance(node, ast.Assign) and assigns:
            if any(isinstance(t, ast.Name) and t.id in assigns for t in node.targets):
                t = text_of(node.value)
                if t:
                    for off, para in _paragraphs(t):
                        out.append((node.value.lineno + off, para))
        elif isinstance(node, ast.Dict) and dict_keys:
            for k, v in zip(node.keys, node.values):
                if isinstance(k, ast.Constant) and k.value in dict_keys:
                    t = text_of(v)
                    if t:
                        out.append((v.lineno, t))
        elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and decorators:
            names = set()
            for d in node.decorator_list:
                f = d.func if isinstance(d, ast.Call) else d
                names.add(name_of(f))
            if names & decorators:
                doc = ast.get_docstring(node, clean=True)
                if doc and node.body and isinstance(node.body[0], ast.Expr):
                    base = node.body[0].lineno
                    for off, para in _paragraphs(doc):
                        out.append((base + off, para))
    return out


def _paragraphs(doc):
    out, buf, start = [], [], 0
    for i, line in enumerate(doc.split("\n")):
        s = line.strip()
        if not s:
            if buf:
                out.append((start, " ".join(buf)))
            buf = []
            continue
        if re.match(r"([-*]|\d+[.)])\s", s) and buf:
            out.append((start, " ".join(buf)))
            buf = []
        if not buf:
            start = i
        buf.append(s)
    if buf:
        out.append((start, " ".join(buf)))
    return out


C_STR = re.compile(r'"((?:[^"\\\n]|\\.)*)"')


def cpp_units(src, cfg, path):
    """String literal runs in help output, log errors and exceptions."""
    src = re.sub(r"/\*.*?\*/", _blank, src, flags=re.S)
    src = re.sub(r"(?m)//.*$", "", src)
    src = re.sub(r"(?m)^\s*#\s*include\b[^\n]*", _blank, src)
    every = any(fnmatch.fnmatch(path, g) for g in cfg.get("cpp_all_strings", []))
    ctx = re.compile(cfg.get("cpp_context", r"(logger->(error|warn|critical)|throw\s+[\w:]+|fprintf\s*\(\s*stderr|std::cerr)"))
    out = []
    pos = 0
    while True:
        m = C_STR.search(src, pos)
        if not m:
            break
        parts, end = [m.group(1)], m.end()
        while True:
            n = re.match(r'\s*"((?:[^"\\\n]|\\.)*)"', src[end:])
            if not n:
                break
            parts.append(n.group(1))
            end += n.end()
        pos = end
        if not every:
            back = src[max(0, m.start() - 300):m.start()]
            hits = list(ctx.finditer(back))
            if not hits or ";" in back[hits[-1].end():]:
                continue
        text = "".join(parts)
        text = text.replace("\\n", " ").replace("\\t", " ").replace('\\"', '"')
        out.append((src.count("\n", 0, m.start()) + 1, text))
    return out


TS_SKIP_ATTR = re.compile(r"(className|class|href|src|id|key|data-[\w-]+|style|type|name|role|variant|size|rel|target|as|method|testId|import|from|require|url|path|route|mode|align|side)\s*[=:(]\s*\{?\s*$")


def ts_units(src):
    """String literals and JSX text in TypeScript and TSX."""
    src = re.sub(r"/\*.*?\*/", _blank, src, flags=re.S)
    src = re.sub(r"(?m)(^|[^:'\"`\\])//.*$", r"\1", src)
    out = []
    pat = re.compile(r"'((?:[^'\\\n]|\\.)*)'|\"((?:[^\"\\\n]|\\.)*)\"|`((?:[^`\\]|\\.)*)`|>([^<>{}]*[A-Za-z][^<>{}]*)<")
    for m in pat.finditer(src):
        s = next(g for g in m.groups() if g is not None)
        before = src[max(0, m.start() - 40):m.start()]
        if m.group(4) is None and TS_SKIP_ATTR.search(before):
            continue
        if re.search(r"\b(console\.\w+|logger\.\w+|log\.\w+|new Error)\s*\([^)]{0,30}$", before):
            continue
        if re.search(r"\bimport\b|\bfrom\s*$", before):
            continue
        s = re.sub(r"\$\{[^}]*\}", "X", s)
        if m.group(4) is not None:
            s = re.sub(r"\s+", " ", s)
        out.append((src.count("\n", 0, m.start()) + 1, s))
    return out


# ---------- checks ----------

def looks_like_code(t):
    if re.search(r"(^|\s)(SELECT|INSERT|UPDATE|DELETE)\s", t):
        return True
    if t.count("(") + t.count("[") + t.count("=") > len(t) / 12:
        return True
    return bool(re.search(r"--?\w[\w-]*=|::|->|\bvar\(|\d+px|\\[dws]|\^\(|\$\(", t))


def sentences(t):
    for s in re.split(r"(?<=[.!?])\s+(?=[A-Z0-9\"'(])", t):
        yield s


def check_text(t, cfg, scope):
    """Return a list of (rule, snippet) for one text unit."""
    t = re.sub(r"\{[^{}]*\}|%[-0-9.]*[sdzlufgx]+|\bCODE\b", "X", t)
    if looks_like_code(t):
        return []
    hits = []
    if ";" in t:
        hits.append(("R1", t))
    for s in sentences(t):
        n = len(WORD.findall(s))
        if n > cfg.get("max_sentence_words", 30):
            hits.append(("R2", f"{n} words: {s}"))
    if re.search(r"\b(sorry|oops|whoops|something went wrong|an error occurred|unknown error)\b", t, re.I):
        hits.append(("R4", t))
    if re.search(r"\b(may|might|soon|shortly|a while|eventually)\b", re.sub(r"\bMay\s+\d", "", t)):
        hits.append(("R7", t))
    elif re.search(r"\b(two|three|four|five|six|seven|eight|nine|ten|twelve|twenty|thirty|sixty)\s+(seconds?|minutes?|hours?|days?|weeks?|months?|times|threads?|files?|steps?|results?|GB|MB|TB)\b", t, re.I):
        hits.append(("R7", t))
    for term in cfg.get("terms", []):
        for bad in term["not"]:
            if re.search(r"(?<![\w-])" + re.escape(bad) + r"(?![\w-])", t, 0 if term.get("case") else re.I):
                hits.append(("R8", f"'{bad}' -> '{term['use']}': {t}"))
    if scope.get("jargon") and PROSE.search(t):  # a lone identifier or a code line is not prose
        for j in cfg.get("jargon", []):
            pat = j[3:] if j.startswith("re:") else r"(?<![\w-])" + re.escape(j) + r"(?![\w-])"
            if re.search(pat, t):
                hits.append(("R5", f"'{j}': {t}"))
    if scope.get("no_we") and re.search(r"\b(we|We|our|Our|us)\b", t):
        hits.append(("R9", t))
    if re.search(r"\b(is|are|was|were|be|been|being)\s+(\w+ly\s+)?(\w+ed|shown|done|made|given|taken|seen|known|built|kept|sent|set|run|written|found|lost|held|read)\b", t):
        hits.append(("R6", t))
    return hits


# ---------- driver ----------

def repo_files(root):
    res = subprocess.run(["git", "-C", str(root), "ls-files"], capture_output=True, text=True, check=True)
    return res.stdout.split("\n")


def scan(root, cfg, only=None):
    files = repo_files(root)
    results = {}
    for scope in cfg["scopes"]:
        for f in files:
            if not f or not any(fnmatch.fnmatch(f, g) for g in scope["include"]):
                continue
            if any(fnmatch.fnmatch(f, g) for g in scope.get("exclude", []) + cfg.get("exclude", [])):
                continue
            if only and f not in only:
                continue
            p = root / f
            try:
                src = p.read_text(errors="replace")
            except (FileNotFoundError, IsADirectoryError):
                continue
            ext = p.suffix
            if ext == ".md":
                units = markdown_units(src)
            elif ext == ".rst":
                units = markdown_units(src, rst=True)
            elif ext == ".py":
                units = python_units(src, cfg)
            elif ext in (".cpp", ".hpp", ".cc", ".h", ".cu", ".cuh"):
                units = cpp_units(src, cfg, f)
            elif ext in (".ts", ".tsx", ".js", ".jsx"):
                units = ts_units(src)
            else:
                continue
            lines = src.split("\n")
            for ln, text in units:
                ctx = " ".join(lines[max(0, ln - 2):ln])
                allowed = set(re.findall(r"R\d+", " ".join(ALLOW.findall(ctx))))
                for rule, snip in check_text(text, cfg, scope):
                    if rule in allowed:
                        continue
                    results.setdefault(f, []).append((ln, rule, re.sub(r"\s+", " ", snip).strip()))
    return results


def counts(results):
    out = {}
    for f, hits in results.items():
        c = {}
        for _, rule, _ in hits:
            if rule in RATCHET:
                c[rule] = c.get(rule, 0) + 1
        if c:
            out[f] = dict(sorted(c.items()))
    return dict(sorted(out.items()))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="*")
    here = pathlib.Path(__file__).resolve().parent
    ap.add_argument("--root", default=str(here.parent))
    ap.add_argument("--config", default=str(here / "text-rules.json"))
    ap.add_argument("--baseline", default=str(here / "text-rules-baseline.json"))
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--update-baseline", action="store_true")
    a = ap.parse_args()
    root = pathlib.Path(a.root).resolve()
    cfg = json.loads(pathlib.Path(a.config).read_text())
    results = scan(root, cfg, set(a.files) or None)
    now = counts(results)

    if a.report:
        totals = {}
        for f in sorted(results):
            for ln, rule, snip in sorted(results[f]):
                print(f"{f}:{ln}: {rule}: {snip[:200]}")
                totals[rule] = totals.get(rule, 0) + 1
        print("TOTALS", json.dumps(dict(sorted(totals.items()))), "files", len(results), file=sys.stderr)
        return 0
    base_path = pathlib.Path(a.baseline)
    base = json.loads(base_path.read_text()) if base_path.exists() else {}
    if a.update_baseline:
        increases = [f"{f}: {rule} has {n} hits, baseline {base.get(f, {}).get(rule, 0)}"
                     for f, rules in now.items() for rule, n in rules.items()
                     if n > base.get(f, {}).get(rule, 0)]
        if increases:
            print("Cannot raise the baseline. Fix these hits or document an exception:")
            print("\n".join(increases))
            return 1
        base_path.write_text(json.dumps(now, indent=1) + "\n")
        print(f"wrote {base_path.name}: {sum(sum(c.values()) for c in now.values())} hits in {len(now)} files")
        return 0
    problems = []
    for f in sorted(set(now) | set(base)):
        for rule in RATCHET:
            n, b = now.get(f, {}).get(rule, 0), base.get(f, {}).get(rule, 0)
            if n > b:
                new = [h for h in results.get(f, []) if h[1] == rule]
                problems.append(f"{f}: {rule} has {n} hits, baseline {b}. The new hit is one of these:")
                problems += [f"  {f}:{ln}: {snip[:160]}" for ln, _, snip in new]
            elif n < b and (f in now or (root / f).exists()):
                problems.append(f"{f}: {rule} has {n} hits, baseline {b}. Run --update-baseline to lower it.")
    if problems:
        print("\n".join(problems))
        print(f"\nRules: see {cfg.get('rules_doc', 'AGENTS.md')}. A deliberate exception takes "
              "`text-rules-allow: <rule>` in a comment on that line or the line before.")
        return 1
    print(f"text rules: OK ({sum(sum(c.values()) for c in now.values())} known hits, none new)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
