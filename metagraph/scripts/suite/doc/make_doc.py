#!/usr/bin/env python3
"""The suite results document: one PDF per run from a recorded run directory (DESIGN-DOC §9).

    make_doc.py RUN_DIR --out DIR [--chapters R,T,...] [--sample] [--no-pdf] [--recorded]
                [--units FILE] [--prose-dir DIR] [--checks-dir DIR] [--by-design PATH]
                [--policy FILE] [--also-policy FILE ...] [--declared FILE] [--identity FILE] [--xelatex PATH]
                [--scaffold-prose]

Reads the run directory (manifest, the re-evaluated executions and compare.json when present, else the recorded
ones, run.json, reevaluated.json, samples, gate.json, arms/mcp/run.json), the unit metadata (the manifest's `doc`
blocks, else the mapping `--units`, default doc/units/suite_v<version>.json), the prose data (`--prose-dir`,
default doc/prose/: <chapter>.json, overview.json, glossary.json, checks.json, closing.json, cases.json) and
writes DIR/doc.tex with one included file per chapter, the plots' data files (plot_*.csv), DIR/doc.sty, then
builds DIR/doc.pdf with xelatex (twice; a third time when the .aux changed). Python stdlib only; offline;
deterministic (no clock in the output; the README's commit is read from the checkout).

--chapters picks the chapters of Part II (the appendices and the index follow); --sample is --chapters R with the
cover note "Chapter R sample; whole-run summary shown for context";
--recorded reads the recorded set instead of the re-evaluated one; --identity names the deployed service's identity
the run did not record (a JSON file: run, service_commit, railway_deployment_id, source, and host_release,
mcp_runner_commit, manifest checked against the run), printed in I.4 and appendix F as taken from its source and
not recorded by the runner; --also-policy (repeatable) judges a further gate policy offline and prints its verdict in
I.4 beside the deploy gate's (--policy, else the first policy of gate/ pinning the manifest), which it does not
replace; --scaffold-prose writes the placeholder
entries the prose files lack (existing text untouched) and exits.

Exit 0 built (or --no-pdf), 1 LaTeX failed (the log's first error echoed), 2 usage / missing inputs (a scenario
without unit metadata, a family or unit without prose, a missing prose file: each named), 4 a comparison to compute
(no compare.json) with a by-design table compare.py refuses (as compare.py and report.py).
"""

import argparse
import json
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
if os.path.dirname(HERE) not in sys.path:
    sys.path.insert(0, os.path.dirname(HERE))

import reader  # noqa: E402
import sections  # noqa: E402
import units as units_mod  # noqa: E402
from latex import esc  # noqa: E402

DEFAULT_XELATEX = "/Library/TeX/texbin/xelatex"
TEMPLATES = os.path.join(HERE, "templates")
DEFAULT_PROSE = os.path.join(HERE, "prose")


class DocError(Exception):
    """Usage or missing inputs (exit 2)."""


def read_json(path):
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)


def load_prose(prose_dir, units, chapters):
    """({chapter: prose}, overview, glossary, checks prose, closing prose, cases prose); a missing file or key for a
    built chapter raises DocError naming it. The shared files overview.json, glossary.json, checks.json (the check
    catalogue), closing.json (the closing section's templates) and cases.json (the short reasons of the by-design
    cases) are optional."""
    prose = {}
    missing = []
    for ch in chapters:
        path = os.path.join(prose_dir, "%s.json" % ch)
        if not os.path.exists(path):
            missing.append("%s (no file %s)" % (ch, path))
            continue
        data = read_json(path)
        prose[ch] = data
        if not data.get("intro"):
            missing.append("%s: intro" % ch)
        fams = data.get("families") or {}
        us = data.get("units") or {}
        for fid in units.chapter_families(ch):
            entry = fams.get(fid) or {}
            if not entry.get("prose"):
                missing.append("%s: families.%s.prose" % (ch, fid))
            for uid in units.families[fid]["units"]:
                if not us.get(uid):
                    missing.append("%s: units.%s" % (ch, uid))
    # the other chapters' prose files feed the check catalogue when present
    for ch in units.chapters:
        if ch in prose:
            continue
        path = os.path.join(prose_dir, "%s.json" % ch)
        if os.path.exists(path):
            prose[ch] = read_json(path)
    if missing:
        raise DocError("prose missing (%d): %s" % (len(missing), "; ".join(missing[:30])
                                                   + (" ..." if len(missing) > 30 else "")))
    shared = {}
    for name in ("overview", "glossary", "checks", "closing", "cases"):
        path = os.path.join(prose_dir, "%s.json" % name)
        shared[name] = read_json(path) if os.path.exists(path) else {}
    return prose, shared["overview"], shared["glossary"], shared["checks"], shared["closing"], shared["cases"]


def render(template, fields):
    out = template
    for key, value in fields.items():
        out = out.replace("{{%s}}" % key, value)
    return out


def build_tex(run, units, chapters, prose, overview, glossary, checks_prose=None, closing_prose=None,
              cases_prose=None, sample=False):
    """(main .tex text, {file name: text}, the context) for the run: the chapter files and the plots' data files
    are written beside the main .tex."""
    ctx = sections.Context(run, units, prose, overview, glossary, chapters, checks_prose, closing_prose,
                           cases_prose, sample)
    with open(os.path.join(TEMPLATES, "main.tex.in"), "r", encoding="utf-8") as fh:
        template = fh.read()
    # the chapters first: Part I's annotated matrix and the decision summary refer to the footnotes they assign
    chapter_files = {}
    part_two = "\\part*{Part II\\qd Results by area}\\addcontentsline{toc}{part}{Part II\\qd Results by area}\n"
    for i, ch in enumerate(chapters):
        name = "chapter_%s.tex" % ch
        chapter_files[name] = sections.chapter(ctx, ch)
        part_two += ("\\clearpage\n" if i else "") + "\\input{%s}\n" % name[:-4]
    part_one = sections.part_one(ctx)
    closing = sections.closing(ctx)
    appendices = sections.appendices(ctx)
    chapter_files.update(ctx.data_files)
    started = None
    for info in (run.run.get("arms") or {}).values():
        s = (info or {}).get("started_utc")
        if s and (started is None or s < started):
            started = s
    fields = {"run_id": esc(run.run_id()), "title": "MetaGraph suite results\\\\[2pt]{\\large\\texttt{%s}}" % esc(
        run.run_id()), "author": "generated by make\\_doc.py from the run directory%s" % (
        " (the re-evaluated set)" if run.reevaluated else " (the recorded set)"),
              "date": esc((started or "")[:10]) or "–", "verdict_line": sections.verdict_line(ctx),
              "gate_line": sections.gate_line(ctx), "sample_note": sections.sample_note(ctx),
              "front_matter": sections.front_matter(ctx), "part_one": part_one, "part_two": part_two,
              "closing": closing, "appendices": appendices}
    return render(template, fields), chapter_files, ctx


def write_outputs(out_dir, main_tex, chapter_files, name="doc"):
    os.makedirs(out_dir, exist_ok=True)
    main_path = os.path.join(out_dir, "%s.tex" % name)
    with open(main_path, "w", encoding="utf-8") as fh:
        fh.write(main_tex)
    for fname, text in chapter_files.items():
        with open(os.path.join(out_dir, fname), "w", encoding="utf-8") as fh:
            fh.write(text)
    shutil.copyfile(os.path.join(TEMPLATES, "doc.sty"), os.path.join(out_dir, "doc.sty"))
    return main_path


def find_xelatex(given=None):
    for candidate in ([given] if given else []) + [DEFAULT_XELATEX, shutil.which("xelatex")]:
        if candidate and os.path.exists(candidate):
            return candidate
    return None


def _read(path):
    try:
        with open(path, "rb") as fh:
            return fh.read()
    except OSError:
        return b""


def build_pdf(out_dir, main_path, xelatex):
    """(ok, message): xelatex twice, a third time when the .aux changed between the passes."""
    name = os.path.splitext(os.path.basename(main_path))[0]
    aux = os.path.join(out_dir, name + ".aux")
    log = os.path.join(out_dir, name + ".log")
    env = dict(os.environ)
    env.pop("TEXINPUTS", None)
    cmd = [xelatex, "-interaction=nonstopmode", "-halt-on-error", "-no-shell-escape", os.path.basename(main_path)]
    previous = None
    for i in range(3):
        if i == 2 and previous == _read(aux):
            break
        previous = _read(aux)
        proc = subprocess.run(cmd, cwd=out_dir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=env)
        if proc.returncode != 0:
            text = _read(log).decode("utf-8", "replace").splitlines()
            lines = []
            for j, l in enumerate(text):
                if l.startswith("!"):
                    lines = text[j:j + 12]
                    break
            return False, "xelatex failed (pass %d): %s" % (i + 1, "\n".join(lines) or proc.stdout.decode(
                "utf-8", "replace")[-1500:])
    pdf = os.path.join(out_dir, name + ".pdf")
    return os.path.exists(pdf), pdf


def scaffold_prose(prose_dir, units, run):
    """Write the placeholder entries every chapter's prose file lacks; existing text untouched."""
    os.makedirs(prose_dir, exist_ok=True)
    checks_path = os.path.join(prose_dir, "checks.json")
    catalogue = (read_json(checks_path) if os.path.exists(checks_path) else {}).get("checks") or {}
    kinds_without_doc = {}
    for sid in run.order:
        ch = (units.by_scenario.get(sid) or {}).get("chapter")
        for e in run.scenario(sid).get("expect") or []:
            name = e.get("check") if isinstance(e, dict) else None
            if name and not run.check_doc(name) and not (catalogue.get(name) or {}).get("asserts"):
                kinds_without_doc.setdefault(name, ch)
    written = []
    for ch in units.chapters:
        path = os.path.join(prose_dir, "%s.json" % ch)
        existing = read_json(path) if os.path.exists(path) else {}
        kinds = sorted(k for k, c in kinds_without_doc.items() if c == ch)
        data = units_mod.prose_scaffold(units, existing, ch, kinds)
        if data != existing:
            with open(path, "w", encoding="utf-8") as fh:
                json.dump(data, fh, indent=1, ensure_ascii=False)
                fh.write("\n")
            written.append(path)
    return written


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir")
    ap.add_argument("--out", required=True)
    ap.add_argument("--chapters", default=None, help="comma-separated chapter letters (default: all)")
    ap.add_argument("--sample", action="store_true", help="Part I, chapter R, the closing skeleton, R's appendices")
    ap.add_argument("--no-pdf", action="store_true")
    ap.add_argument("--recorded", action="store_true", help="the recorded set instead of the re-evaluated one")
    ap.add_argument("--units", default=None, help="the unit mapping (default doc/units/suite_v<version>.json)")
    ap.add_argument("--prose-dir", default=DEFAULT_PROSE)
    ap.add_argument("--checks-dir", default=None)
    ap.add_argument("--by-design", default=None)
    ap.add_argument("--policy", default=None)
    ap.add_argument("--also-policy", action="append", default=[], help="a further gate policy judged offline and "
                    "shown in I.4 beside the deploy gate's (repeatable)")
    ap.add_argument("--declared", default=None)
    ap.add_argument("--identity", default=None, help="the deployed service's identity (a JSON file the run did not "
                    "record, e.g. from the deployment log)")
    ap.add_argument("--xelatex", default=None)
    ap.add_argument("--name", default="doc", help="the base name of the .tex and .pdf (default doc)")
    ap.add_argument("--scaffold-prose", action="store_true", help="write the missing prose placeholders and exit")
    args = ap.parse_args(argv)
    try:
        run = reader.Run(args.run_dir, recorded=args.recorded, by_design=args.by_design, checks_dir=args.checks_dir,
                         declared=args.declared, policy=args.policy, identity=args.identity,
                         also_policies=args.also_policy)
        mapping = None
        needs_mapping = any(not isinstance(s.get("doc"), dict) for s in run.scenarios.values())
        if needs_mapping:
            path = args.units or units_mod.default_mapping_path(run.manifest)
            if not os.path.exists(path):
                raise DocError("no unit mapping at %s and the manifest has no doc blocks" % path)
            mapping = units_mod.load_mapping(path)
        units = units_mod.resolve(run.manifest, mapping)
        problems = units_mod.validate(units, run.manifest)
        if problems:
            raise DocError("unit metadata (%d problems): %s" % (len(problems), "; ".join(problems[:20])))
        if args.scaffold_prose:
            for p in scaffold_prose(args.prose_dir, units, run):
                print("wrote %s" % p)
            return 0
        chapters = ["R"] if args.sample else ([c.strip() for c in args.chapters.split(",") if c.strip()]
                                             if args.chapters else list(units.chapters))
        unknown = [c for c in chapters if c not in units.chapters]
        if unknown:
            raise DocError("unknown chapter(s): %s (the run has %s)" % (", ".join(unknown), ", ".join(units.chapters)))
        prose, overview, glossary, checks_prose, closing_prose, cases_prose = load_prose(args.prose_dir, units,
                                                                                           chapters)
        main_tex, chapter_files, ctx = build_tex(run, units, chapters, prose, overview, glossary, checks_prose,
                                                 closing_prose, cases_prose, sample=args.sample)
    except reader.TableRefused as e:
        print("make_doc: refusing: %s" % e, file=sys.stderr)
        return 4
    except (reader.RunError, units_mod.UnitsError, DocError) as e:
        print("make_doc: %s" % e, file=sys.stderr)
        return 2
    main_path = write_outputs(args.out, main_tex, chapter_files, args.name)
    print("wrote %s (%d chapter file(s))" % (main_path, len(chapter_files)))
    unresolved = sorted(set(ctx.unresolved))
    if unresolved:
        # a placeholder or condition the run's values lack stays visible in the document; say so here too
        print("make_doc: WARNING: %d prose placeholder(s) or condition(s) the run's values lack (printed as "
              "written): %s" % (len(unresolved), ", ".join(unresolved[:40]) + (" ..." if len(unresolved) > 40 else "")),
              file=sys.stderr)
    if args.no_pdf:
        return 0
    xelatex = find_xelatex(args.xelatex)
    if not xelatex:
        print("make_doc: xelatex not found (looked at %s and PATH); the .tex is in %s" % (DEFAULT_XELATEX, args.out),
              file=sys.stderr)
        return 1
    ok, message = build_pdf(args.out, main_path, xelatex)
    if not ok:
        print("make_doc: %s" % message, file=sys.stderr)
        return 1
    print(message)
    return 0


if __name__ == "__main__":
    sys.exit(main())
