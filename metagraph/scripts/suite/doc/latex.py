"""LaTeX emission for the suite results document (DESIGN-DOC §4, §9.3): escaping, the small typesetting helpers, the
five unit-table emitters (series, replicates, refusals, list, pairs) and the detail-block box.

Every string that comes from a run (ids, titles, sequences, error texts, check details) passes through `esc`: the
ten LaTeX specials are escaped, control characters dropped, sequences of 41+ bases cut to `ACGT...(123 bp)` as
report.ascii_safe does, and the symbols the text fonts may lack are set in the glyph font (DejaVu Sans). Nothing
here reads a run; the emitters take ready strings and return LaTeX.
"""

import re

TEXT_WIDTH_MM = 174.0  # A4, 18 mm margins
TABCOLSEP_MM = 1.06  # 3 pt, each side of a column

_SPECIALS = {"\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "$": r"\$", "#": r"\#", "_": r"\_", "{": r"\{",
             "}": r"\}", "~": r"\textasciitilde{}", "^": r"\textasciicircum{}"}
_SPECIAL_RE = re.compile(r"[\\&%$#_{}~^]")
_SEQ_RE = re.compile(r"[ACGTNacgtnRYSWKMBDHVU]{41,}")
_CTRL_RE = re.compile(r"[\x00-\x08\x0b-\x1f\x7f]")
# code points the text fonts (Libertinus Serif, Fira Sans, Fira Mono) may lack: set in the glyph font
_SYMBOLS = set("≠≤≥⊆⊇⊑⊒∈∉∅≈⊘●✖▲◐○◆→←⇒⊂⊃∩∪")


def cut_sequences(text):
    """Sequences of 41 or more bases cut to their first 16 and their length, as report.ascii_safe does."""
    return _SEQ_RE.sub(lambda m: "%s...(%d bp)" % (m.group(0)[:16], len(m.group(0))), text)


def cut_text(text, limit):
    """`text` cut to at most `limit` characters with an ellipsis, at a word boundary when one lies in the last
    two fifths of the kept part (never mid-word otherwise); unchanged when it fits."""
    if limit is None or len(text) <= limit:
        return text
    head = text[:max(0, limit - 1)]
    space = max(head.rfind(" "), head.rfind("\u00a0"))
    if space >= int(limit * 0.6):
        head = head[:space]
    return head.rstrip(" ,;:(\u00a0") + "…"


def esc(text, limit=None):
    """`text` safe inside LaTeX text: specials escaped, control characters and line breaks turned into spaces, long
    sequences cut, symbols the fonts may lack set in the glyph font, cut at `limit` characters with an ellipsis."""
    if text is None:
        return ""
    t = str(text).replace("\r\n", " ").replace("\n", " ").replace("\r", " ")
    t = _CTRL_RE.sub("", t)
    t = cut_sequences(t)
    t = cut_text(t, limit)
    out = []
    for ch in t:
        if ch in _SPECIALS:
            out.append(_SPECIALS[ch])
        elif ch in _SYMBOLS:
            out.append("{\\mdseries\\upshape\\glyphfont %s}" % ch)
        else:
            out.append(ch)
    return "".join(out)


_BREAK_CHARS = "/.-_:"
_LONG_RUN_RE = re.compile(r"[A-Za-z0-9]{14,}")


def _breakable(text, limit=None):
    """`text` escaped, with a break opportunity after each of / . - _ : and inside long alphanumeric runs (hashes),
    so that ids, URLs and case ids wrap inside narrow monospace cells."""
    raw = str(text) if text is not None else ""
    raw = _CTRL_RE.sub("", raw.replace("\r\n", " ").replace("\n", " ").replace("\r", " "))
    raw = cut_sequences(raw)
    raw = cut_text(raw, limit)
    raw = _LONG_RUN_RE.sub(lambda m: "\x00".join(m.group(0)[i:i + 10] for i in range(0, len(m.group(0)), 10)), raw)
    out = []
    for ch in raw:
        if ch == "\x00":
            out.append("\\allowbreak{}")
            continue
        out.append(esc(ch))
        if ch in _BREAK_CHARS:
            out.append("\\allowbreak{}")
    return "".join(out)


_TOKEN_RE = re.compile(r"\S{18,}")


def esc_wrap(text, limit=None):
    """`text` escaped as `esc` does, with a break opportunity after . _ / - inside any token of 18+ characters
    (case ids and field paths quoted in free text), for the narrow columns of the appendix tables."""
    if text is None:
        return ""
    t = str(text).replace("\r\n", " ").replace("\n", " ").replace("\r", " ")
    t = _CTRL_RE.sub("", t)
    t = cut_sequences(t)
    t = cut_text(t, limit)
    out = []
    pos = 0
    for m in _TOKEN_RE.finditer(t):
        out.append(esc(t[pos:m.start()]))
        token = m.group(0)
        piece = ""
        for ch in token:
            piece += esc(ch)
            if ch in "._/-":
                piece += "\\allowbreak{}"
        out.append(piece)
        pos = m.end()
    out.append(esc(t[pos:]))
    return "".join(out)


def tt(text, limit=None):
    """Monospace text (an id, a code, a quoted error text), breakable after punctuation."""
    return "\\texttt{%s}" % _breakable(text, limit)


def title_tt(text):
    """Monospace text in a heading, with a plain PDF-string twin for the bookmarks."""
    return "\\texorpdfstring{\\texttt{%s}}{%s}" % (esc(text), esc(text))


def label_key(text):
    """A label or reference key: ASCII letters, digits, ':' and '-' only."""
    return re.sub(r"[^A-Za-z0-9:-]", "-", str(text))


def breakable_id(text):
    """An id in monospace, breakable after its hyphens and underscores (longtable columns are narrow)."""
    return tt(text)


def sid(scenario_id, unit_id=None, detail=False):
    """A test id, linked to its unit's table (and, with `detail`, to its detail block)."""
    body = breakable_id(scenario_id)
    if detail:
        return "\\hyperref[detail:%s]{%s}" % (label_key(scenario_id), body)
    if unit_id:
        return "\\hyperref[unit:%s]{%s}" % (label_key(unit_id), body)
    return body


_UNIT_RE = re.compile(r"(?<![A-Za-z\d.,])(\d[\d,]*(?:\.\d+)?) (ms|s|bp|kb|MiB|GiB|\\%)(?![A-Za-z])")


def join_units(text):
    """An escaped prose string with a non-breaking space between a number and its unit (`25 s`, `2,000 bp`,
    `7762 \\%`), so that a filled figure never breaks from its unit at a line's end (re-verification of round 3,
    low item 8)."""
    return _UNIT_RE.sub(lambda m: "%s~%s" % (m.group(1), m.group(2)), text)


def inline(text):
    """Prose with the tiny inline markup of the prose files: `code` and *emphasis*; everything else escaped, a
    number joined to its unit by a non-breaking space (join_units)."""
    parts = re.split(r"(`[^`]*`|\*[^*]+\*)", str(text or ""))
    out = []
    for p in parts:
        if len(p) >= 2 and p.startswith("`") and p.endswith("`"):
            code = p[1:-1]
            out.append(tt(code) if (len(code) > 24 or "/" in code) else "\\texttt{%s}" % esc(code))
        elif len(p) >= 2 and p.startswith("*") and p.endswith("*"):
            out.append("\\emph{%s}" % join_units(esc(p[1:-1])))
        else:
            out.append(join_units(esc(p)))
    return "".join(out)


def paragraphs(items):
    """Prose paragraphs (a list of strings, or one string) as LaTeX paragraphs."""
    if isinstance(items, str):
        items = [items]
    return "\n\n".join(inline(p) for p in (items or []) if p) + ("\n\n" if items else "")


def placeholder(text):
    """A prose placeholder, visibly marked."""
    return "{\\color{black!55}\\itshape %s}" % esc(text)


def col_p(width_mm, align="raggedright"):
    return ">{\\%s\\arraybackslash}p{%.2fmm}" % (align, width_mm)


def fit(widths, extra_cols=0, extra_mm=15.0):
    """Column widths scaled so that a table of these p-columns plus `extra_cols` natural-width columns (numbers
    under a header of up to ten letters, `extra_mm` each) and @{} at both ends fits the text width with the 3 pt
    column separation."""
    n = len(widths) + extra_cols
    available = TEXT_WIDTH_MM - 2 * TABCOLSEP_MM * (n - 1) - extra_mm * extra_cols
    total = float(sum(widths)) or 1.0
    return [w * available / total for w in widths]


def spec_p(widths, extra_cols=0, aligns=None, extra_spec=""):
    """A column spec of p-columns that fit beside `extra_spec` (the natural-width columns, e.g. 'r r')."""
    aligns = aligns or ["raggedright"] * len(widths)
    cols = " ".join(col_p(w, a) for w, a in zip(fit(widths, extra_cols), aligns))
    return "@{}%s%s@{}" % (cols, (" " + extra_spec) if extra_spec else "")


def cell_lines(line1, line2=None, line3=None):
    """A table cell of up to three lines: the second and third in small grey type."""
    out = [line1 or "~"]
    for extra in (line2, line3):
        if extra:
            out.append("{\\scriptsize\\color{black!60}%s}" % extra)
    return "\\newline ".join(out)


def _ids_row(ncols, text):
    return "\\multicolumn{%d}{@{}p{%.0fmm}@{}}{\\scriptsize %s}\\\\" % (ncols, TEXT_WIDTH_MM, text)


def unit_header(unit_id, title, unit_type, property_tex, needspace=20):
    """The line above a unit's table: id, title, layout; then the property line; then the label."""
    return ("\\Needspace{%d\\baselineskip}\n\\unithead{%s}{%s}{%s}\n\\propline{%s}\\label{unit:%s}\n" % (
        needspace, esc(unit_id), title, esc(unit_type), property_tex, label_key(unit_id)))


def series_table(columns, rows, ids_text, legend=None):
    """A series matrix: arms as rows, variants as columns. `columns`: the variant labels (LaTeX); `rows`:
    [(arm label, [cell LaTeX, ...])]; `ids_text`: the id line (LaTeX); `legend`: the timing legend line."""
    n = max(1, len(columns))
    label_w = 20.0
    w = (TEXT_WIDTH_MM - label_w - 2 * TABCOLSEP_MM * (n + 1)) / n
    lines = ["\\begingroup\\sffamily\\footnotesize\\setlength\\tabcolsep{3pt}",
             "\\noindent\\begin{tabular}{@{}p{%.1fmm}*{%d}{%s}@{}}" % (label_w, n, col_p(w)),
             "\\toprule",
             " & " + " & ".join(columns) + "\\\\",
             "\\midrule"]
    for arm, cells in rows:
        lines.append("%s & %s\\\\" % (arm, " & ".join(c or "~" for c in cells)))
    lines += ["\\bottomrule", _ids_row(n + 1, ids_text)]
    if legend:
        lines.append(_ids_row(n + 1, "\\color{black!60}" + legend))
    lines += ["\\end{tabular}", "\\endgroup", ""]
    return "\n".join(lines)


def replicates_table(unit_label, rows, ids_text, legend=None):
    """Rows: [(group label, n, cell dl, cell dr, cell mcp)]."""
    spec = "@{}%s r %s %s %s@{}" % (col_p(26), col_p(43), col_p(43), col_p(43))
    header = "%s & n & direct-local & direct-remote & mcp" % unit_label
    return _longtable_lines(spec, header, _rows_with_notes([[g, str(n), a, b, c] for g, n, a, b, c in rows], 5),
                            ids_text, 5, legend)


def _rows_with_notes(rows, ncols):
    """Rows as LaTeX lines; a row given as (cells, note) gets the note in a full-width line under it (the owner's
    review, point 5: a short explanation under the row rather than a smaller font)."""
    out = []
    for r in rows:
        note = None
        if isinstance(r, tuple) and len(r) == 2 and isinstance(r[0], list):
            r, note = r
        out.append(" & ".join(c or "~" for c in r) + "\\\\")
        if note:
            out.append("\\multicolumn{%d}{@{}p{%.0fmm}@{}}{\\hspace*{4mm}%s}\\\\" % (ncols, TEXT_WIDTH_MM, note))
    return out


def _longtable_lines(spec, header, body_lines, ids_text, ncols, legend=None):
    lines = ["\\begingroup\\sffamily\\footnotesize\\setlength\\tabcolsep{3pt}",
             "\\setlength\\LTleft{0pt}\\setlength\\LTright{0pt}",
             "\\begin{longtable}{%s}" % spec,
             "\\toprule", header + "\\\\", "\\midrule", "\\endhead",
             "\\midrule", "\\multicolumn{%d}{r}{\\scriptsize\\itshape continued}\\\\" % ncols, "\\endfoot",
             "\\bottomrule"]
    if ids_text:
        lines.append(_ids_row(ncols, ids_text))
    if legend:
        lines.append(_ids_row(ncols, "\\color{black!60}" + legend))
    lines.append("\\endlastfoot")
    lines += body_lines
    lines += ["\\end{longtable}", "\\endgroup", ""]
    return "\n".join(lines)


def refusals_table(rows, ids_text=None, legend=None):
    """Rows: ([id, input, host expectation, glyph dl, glyph dr, mcp cell], note or None). The mcp column holds the
    recorded status, the service's code and the class of the outcome in words."""
    spec = "@{}%s %s %s %s %s %s@{}" % (col_p(22), col_p(28), col_p(56), col_p(5, "centering"),
                                          col_p(5, "centering"), col_p(44))
    header = "id & input & host expectation & dl & dr & mcp"
    return _longtable_lines(spec, header, _rows_with_notes(rows, 6), ids_text, 6, legend)


def list_table(rows, ids_text=None, legend=None):
    """Rows: ([id, input, host expectation, cell dl, cell dr, cell mcp], note or None)."""
    spec = "@{}%s %s %s %s %s %s@{}" % (col_p(21), col_p(26), col_p(40), col_p(21), col_p(21), col_p(34))
    header = "id & input & host expectation & direct-local & direct-remote & mcp"
    return _longtable_lines(spec, header, _rows_with_notes(rows, 6), ids_text, 6, legend)


def pairs_table(rows, ids_text=None, second="test", legend=None):
    """Rows: [(id, test or source, cell dl, cell dr, cell mcp)]."""
    spec = "@{}%s %s %s %s %s@{}" % (col_p(30), col_p(48), col_p(28), col_p(28), col_p(28))
    header = "id & %s & direct-local & direct-remote & mcp" % second
    return _longtable_lines(spec, header, _rows_with_notes(rows, 5), ids_text, 5, legend)


def plain_table(spec, header, rows, caption=None, small=True, needspace=None):
    """A plain longtable for Part I, the closing section and the appendices (cells already LaTeX); a row may be
    (cells, note) with a full-width note line under it. `needspace` lines are asked for before a table that
    should not start at a page's bottom (default 6 with a caption)."""
    ncols = len(header)
    lines = ["\\begingroup\\sffamily%s\\setlength\\tabcolsep{3pt}" % ("\\footnotesize" if small else "\\small")]
    need = needspace if needspace is not None else (6 if caption else None)
    if need:
        lines.append("\\Needspace{%d\\baselineskip}" % need)
    if caption:
        lines.append("\\noindent{\\itshape %s}\\par\\nopagebreak\\vspace{-4pt}" % caption)
    lines += ["\\setlength\\LTleft{0pt}\\setlength\\LTright{0pt}",
              "\\begin{longtable}{%s}" % spec]
    lines += ["\\toprule", " & ".join(header) + "\\\\", "\\midrule", "\\endhead",
              "\\midrule", "\\multicolumn{%d}{r}{\\scriptsize\\itshape continued}\\\\" % ncols, "\\endfoot",
              "\\bottomrule", "\\endlastfoot"]
    if not rows:
        lines.append("\\multicolumn{%d}{l}{\\itshape none}\\\\" % ncols)
    lines += _rows_with_notes(rows, ncols)
    lines += ["\\end{longtable}", "\\endgroup", ""]
    return "\n".join(lines)


def detail_box(scenario_id, number, title, body, long=False):
    """A test's detail block, labelled detail:<id>: a short block is one unbreakable box (moved whole to the
    next page when it does not fit), a long block breaks across pages and repeats the test id on every
    continuation (pagination: no orphaned continuation of a short block)."""
    head = "\\fn{%d}\\ %s\\quad %s" % (number, breakable_id(scenario_id), title)
    cont = "%s\\ (continued)" % breakable_id(scenario_id)
    env = "detaillong" if long else "detailshort"
    return ("\\phantomsection\\label{detail:%s}\n\\begin{%s}[title={%s},title after break={%s}]\n%s\n\\end{%s}\n" % (
        label_key(scenario_id), env, head, cont, body, env))


def block_table(header, rows):
    """The compact execution table of a detail block: one row per arm and repeat."""
    spec = "@{}l l l l l l@{}"
    lines = ["\\par\\smallskip\\noindent\\begingroup\\sffamily\\footnotesize\\setlength\\tabcolsep{4pt}",
             "\\begin{tabular}{%s}" % spec, "\\toprule", " & ".join(header) + "\\\\", "\\midrule"]
    for r in rows:
        lines.append(" & ".join(c or "~" for c in r) + "\\\\")
    lines += ["\\bottomrule", "\\end{tabular}", "\\endgroup\\par\\smallskip"]
    return "\n".join(lines)


def finding_line(text):
    """The evidence of a finding, set apart from an allowed 503 (an orange rule and the word FINDING)."""
    return "\\findingline{%s}" % text


def bline(key, text):
    """A line of a detail block: a bold key in the margin, hanging text."""
    return "\\bline{%s}{%s}" % (key, text)


def bsub(key, text):
    """An indented sub-line of a detail block (one check, one repeat)."""
    return "\\bsub{%s}{%s}" % (key, text)


def glossary_entry(term, text, see=None):
    """One glossary entry: the term in bold, the text hanging under it, the cross-references in small grey type."""
    out = "\\par\\noindent\\hangindent=3mm\\hangafter=1\\textbf{%s}\\quad %s" % (term, text)
    if see:
        out += " {\\color{black!60}(see %s)}" % see
    return out + "\\par\\vspace{1.5pt}\n"


def catline(key, text):
    """A line of the check catalogue: a bold key in the margin, hanging text, in small type."""
    return "\\catline{%s}{%s}" % (key, text)


def section(title, label=None, star=True, level="section", mark=None, needspace=8):
    """A heading outside the numbering, in the table of contents; `mark` sets the page header. A heading asks for
    `needspace` lines below it, so that it is never left at a page's bottom (pagination)."""
    out = "\\Needspace{%d\\baselineskip}\n\\%s%s{%s}" % (needspace, level, "*" if star else "", title)
    if star:
        out += "\\addcontentsline{toc}{%s}{%s}" % (level, title)
    if label:
        out += "\\label{%s}" % label_key(label)
    if mark:
        out += "\\markright{%s}" % esc(mark)
    return out + "\n"


def tlegend(text):
    """A local legend line under a table (what / and the brackets mean there)."""
    return "\\noindent{\\sffamily\\scriptsize\\color{black!60}%s}\\par\\smallskip\n" % text


PAIRED_PLOT_WIDTH_MM = 62  # the plotting area of a paired scatter (the owner's second review, point 1: 60-70 mm)
PAIRED_PLOT_HEIGHT_MM = 55


def paired_plot(data_file, xcol, ycol, xlabel, ylabel, label):
    """A log-log scatter of paired walls (one point per scenario) with the diagonal y = x, from a data file the
    generator wrote (deterministic; pgfplots reads the columns by name). The plotting area itself is
    PAIRED_PLOT_WIDTH_MM by PAIRED_PLOT_HEIGHT_MM (`scale only axis`: the labels and ticks sit outside it), so that
    two plots side by side each show their distribution at about 62 mm, not a quarter of the text width."""
    return ("\\begin{tikzpicture}\n"
            "\\begin{axis}[scale only axis,width=%dmm,height=%dmm,xmode=log,ymode=log,"
            "xlabel={%s},ylabel={%s},tick label style={font=\\footnotesize},label style={font=\\footnotesize},"
            "grid=major,grid style={black!10},log ticks with fixed point]\n"
            "\\addplot[only marks,mark=*,mark size=0.7pt,color=bdblue,opacity=0.55] table[x=%s,y=%s,col sep=comma]{%s};\n"
            "\\addplot[domain=1:100000,samples=2,color=black!50,dashed] {x};\n"
            "\\end{axis}\\end{tikzpicture}\\label{%s}\n" % (PAIRED_PLOT_WIDTH_MM, PAIRED_PLOT_HEIGHT_MM, xlabel, ylabel,
                                                              xcol, ycol, data_file, label_key(label)))


def figure_block(body, caption):
    """A plot (or two side by side) and its caption kept on one page (a minipage moves as a whole)."""
    return ("\\par\\medskip\\noindent\\begin{minipage}{\\textwidth}\\centering\n%s\n"
            "\\par\\smallskip\\begin{flushleft}{\\small\\itshape %s}\\end{flushleft}\\end{minipage}\\par\\medskip\n" % (
                body, caption))


def phase_boxplot(boxes, label):
    """Box plots of the MCP phases (one box per phase: p10, p25, p50, p75, p90), log scale, from the quantiles
    the generator computed; `boxes`: [(name, p10, p25, p50, p75, p90, n[, n_zero])] in ms. A tick label is two
    lines, the phase and its n (`n 425 of 547` for a phase whose zeros the log scale leaves out: the caption says
    why), centred under its box, so that neighbouring labels never overlap (re-verification N-2)."""
    def _label(b):
        zeros = b[7] if len(b) > 7 else 0
        return "{%s\\\\(n %d%s)}" % (b[0], b[6], (" of %d" % (b[6] + zeros)) if zeros else "")
    lines = ["\\begin{tikzpicture}",
             "\\begin{axis}[width=0.95\\textwidth,height=0.38\\textwidth,ymode=log,boxplot/draw direction=y,"
             "ylabel={ms},tick label style={font=\\footnotesize},label style={font=\\footnotesize},"
             "xticklabel style={align=center,font=\\footnotesize},"
             "xtick={%s},xticklabels={%s},grid=major,grid style={black!10},log ticks with fixed point,"
             "boxplot/box extend=0.5]" % (",".join(str(i + 1) for i in range(len(boxes))),
                                           ",".join(_label(b) for b in boxes))]
    for b in boxes:
        lines.append("\\addplot[boxplot prepared={lower whisker=%.4g,lower quartile=%.4g,median=%.4g,"
                     "upper quartile=%.4g,upper whisker=%.4g},fill=bdblue!15,draw=bdblue] coordinates {};" % (
                         b[1], b[2], b[3], b[4], b[5]))
    lines += ["\\end{axis}\\end{tikzpicture}\\label{%s}" % label_key(label), ""]
    return "\n".join(lines)
