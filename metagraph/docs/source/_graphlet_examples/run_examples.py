#!/usr/bin/env python3
"""Run the code examples of graphlets.rst against the committed fixtures, and check or
rewrite the outputs printed beside them.

    python3 run_examples.py                    # check: each example reproduces its output
    python3 run_examples.py --write            # rewrite the outputs in graphlets.rst
    python3 run_examples.py --write --pending  # also the pending ones, once the fixes land
    python3 run_examples.py --show ID          # print one example's code and actual output

This is not a test of the library (that is api/python/tests); it keeps the page honest.

What counts as an example: a ``.. code-block:: python`` block that follows a reST comment
``.. graphlet-example: <id> [flags]``. Its output is the ``.. code-block:: text`` block
that is the next explicit markup after the python block (prose may sit in between); an
example that prints nothing has no output block. All examples of the page run in order
in ONE namespace, like an interactive session, with the working directory set to the
fixture directory api/python/tests/data/traverse/documents and api/python on sys.path.

Flags:
  pending  the output shown is the intended behaviour of a library fix that is not merged
           yet. A mismatch is reported but does not fail the check, and --write leaves the
           shown output alone unless --pending is given (run it once the fixes are in).
  server   the example talks to a MetaGraph server. The runner starts a stand-in on a free
           local port that answers GET /traverse/capabilities and POST /traverse with the
           committed fixture response whose request has the same seed, and substitutes that
           port for the literal 5555 in the code.

Exit status: 0 when every non-pending example ran and matched (or was rewritten), 1
otherwise. Requires Python >= 3.10 (the library's dataclasses use slots).
"""

import argparse
import contextlib
import copy
import http.server
import io
import json
import os
import re
import sys
import threading
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
SOURCE = os.path.dirname(HERE)
PAGE = os.path.join(SOURCE, 'graphlets.rst')
METAGRAPH = os.path.dirname(os.path.dirname(SOURCE))          # .../metagraph
API = os.path.join(METAGRAPH, 'api', 'python')
FIXTURES = os.path.join(API, 'tests', 'data', 'traverse', 'documents')

MARKER = re.compile(r'^\.\. graphlet-example: ([A-Za-z0-9_.-]+)((?: [a-z]+)*)\s*$')
INDENT = '   '


class Example:
    def __init__(self, ident, flags, code, code_span, out, out_span):
        self.id = ident
        self.flags = flags
        self.code = code
        self.code_span = code_span          # (first body line, end) in the page's lines
        self.out = out                      # the shown output, None without a block
        self.out_span = out_span
        self.actual = None
        self.error = None

    @property
    def pending(self):
        return 'pending' in self.flags


def _block(lines, start):
    """The body of the directive at lines[start]: (first, end, dedented text). The body is
    every following line that is blank or indented by INDENT; trailing blanks are left
    outside it."""
    i = start + 1
    end = i
    while i < len(lines) and (not lines[i].strip() or lines[i].startswith(INDENT)):
        if lines[i].strip():
            end = i + 1
        i += 1
    first = start + 1
    while first < end and not lines[first].strip():
        first += 1
    body = [l[len(INDENT):] if l.startswith(INDENT) else '' for l in lines[first:end]]
    return first, end, '\n'.join(body)


def parse_page(lines):
    examples = []
    i = 0
    while i < len(lines):
        m = MARKER.match(lines[i])
        if not m:
            i += 1
            continue
        ident, flags = m.group(1), set(m.group(2).split())
        j = i + 1
        while j < len(lines) and lines[j].strip() != '.. code-block:: python':
            if MARKER.match(lines[j]):
                raise SystemExit('example %r has no python block' % ident)
            j += 1
        if j == len(lines):
            raise SystemExit('example %r has no python block' % ident)
        cf, ce, code = _block(lines, j)
        k = ce
        while k < len(lines) and not lines[k].startswith('.. ') and not MARKER.match(lines[k]):
            k += 1
        out = out_span = None
        if k < len(lines) and lines[k].strip() == '.. code-block:: text':
            of, oe, out = _block(lines, k)
            out_span = (of, oe)
            k = oe
        examples.append(Example(ident, flags, code, (cf, ce), out, out_span))
        i = k
    ids = [e.id for e in examples]
    dup = {x for x in ids if ids.count(x) > 1}
    if dup:
        raise SystemExit('duplicate example ids: %s' % ', '.join(sorted(dup)))
    return examples


# ------------------------------------------------------------------ the stand-in server

def _fixture_responses():
    by_seq = {}
    for name in sorted(os.listdir(FIXTURES)):
        if not name.endswith('.request.json'):
            continue
        stem = name[:-len('.request.json')]
        with open(os.path.join(FIXTURES, name)) as f:
            req = json.load(f)
        with open(os.path.join(FIXTURES, stem + '.graphlet.json')) as f:
            resp = json.load(f)
        mode = req['strategy'].get('labels', {}).get('mode', 'constrain')
        by_seq.setdefault((req['seeds'][0]['sequence'], mode), resp)
    return by_seq


class _Handler(http.server.BaseHTTPRequestHandler):
    def log_message(self, *args):
        pass

    def _send(self, status, obj):
        data = json.dumps(obj, separators=(',', ':')).encode()
        self.send_response(status)
        self.send_header('Content-Type', 'application/json')
        self.send_header('Content-Length', str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def do_GET(self):
        if self.path.rstrip('/').endswith('/traverse/capabilities'):
            any_resp = next(iter(self.server.responses.values()))
            return self._send(200, any_resp['capabilities'])
        self._send(404, {'error': 'no route %s' % self.path})

    def do_POST(self):
        body = json.loads(self.rfile.read(int(self.headers['Content-Length'])))
        if not self.path.rstrip('/').endswith('/traverse'):
            return self._send(404, {'error': 'no route %s' % self.path})
        seq = body['seeds'][0]['sequence']
        mode = ((body.get('strategy') or {}).get('labels') or {}).get('mode', 'constrain')
        resp = self.server.responses.get((seq, mode))
        if resp is None:
            return self._send(400, {'error': 'the stand-in only replays the committed '
                                             'fixture responses'})
        self._send(200, copy.deepcopy(resp))


@contextlib.contextmanager
def stand_in_server():
    server = http.server.ThreadingHTTPServer(('127.0.0.1', 0), _Handler)
    server.responses = _fixture_responses()
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield server.server_address[1]
    finally:
        server.shutdown()
        server.server_close()


# ------------------------------------------------------------------ running

def _norm(text):
    lines = [l.rstrip() for l in (text or '').split('\n')]
    while lines and not lines[-1]:
        lines.pop()
    return '\n'.join(lines)


def run(examples):
    if API not in sys.path:
        sys.path.insert(0, API)
    ns = {'__name__': '__graphlet_examples__'}
    cwd = os.getcwd()
    os.chdir(FIXTURES)
    try:
        with stand_in_server() as port:
            for ex in examples:
                code = ex.code
                if 'server' in ex.flags:
                    code = code.replace('5555', str(port))
                buf = io.StringIO()
                try:
                    with contextlib.redirect_stdout(buf):
                        exec(compile(code, '<example %s>' % ex.id, 'exec'), ns)
                except Exception:
                    ex.error = traceback.format_exc(limit=4)
                ex.actual = _norm(buf.getvalue())
    finally:
        os.chdir(cwd)


def rewrite(lines, examples, include_pending):
    """The page with the selected examples' output blocks replaced (bottom up, so that
    earlier spans stay valid)."""
    out = list(lines)
    for ex in sorted(examples, key=lambda e: e.code_span[0], reverse=True):
        if ex.error or (ex.pending and not include_pending) or ex.out_span is None:
            continue
        if _norm(ex.out) == ex.actual:
            continue
        body = [INDENT + l if l else '' for l in ex.actual.split('\n')]
        a, b = ex.out_span
        out[a:b] = body
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--write', action='store_true', help='rewrite the output blocks')
    ap.add_argument('--pending', action='store_true',
                    help='with --write: also rewrite the outputs of pending examples')
    ap.add_argument('--show', metavar='ID', help='print one example\'s code and output')
    ap.add_argument('--page', default=PAGE, help='the reST page (default: graphlets.rst)')
    args = ap.parse_args(argv)

    with open(args.page, encoding='utf-8') as f:
        lines = f.read().split('\n')
    examples = parse_page(lines)
    run(examples)

    if args.show:
        ex = next((e for e in examples if e.id == args.show), None)
        if ex is None:
            raise SystemExit('no example %r (known: %s)' % (args.show,
                                                           ', '.join(e.id for e in examples)))
        print(ex.code)
        print('--- output' + (' (pending)' if ex.pending else ''))
        print(ex.actual)
        if ex.error:
            print('--- error\n' + ex.error)
        return 0

    failed = 0              # errors, missing output blocks, unwritten mismatches
    pending = 0
    for ex in examples:
        tag = 'pending ' if ex.pending else ''
        if ex.error:
            print('ERROR   %s%s\n%s' % (tag, ex.id, ex.error))
            failed += not ex.pending
            pending += ex.pending
            continue
        if ex.out_span is None:
            if ex.actual:
                print('MISSING %s%s: prints output but has no text block' % (tag, ex.id))
                failed += not ex.pending
                pending += ex.pending
            else:
                print('ok      %s%s' % (tag, ex.id))
            continue
        if _norm(ex.out) == ex.actual:
            print('ok      %s%s' % (tag, ex.id))
        elif args.write and (not ex.pending or args.pending):
            print('written %s%s' % (tag, ex.id))
        else:
            print('DIFFERS %s%s' % (tag, ex.id))
            if ex.pending:
                print('        (pending: the shown output is the intended behaviour; '
                      'the library prints today:)')
            for l in ex.actual.split('\n'):
                print('        | ' + l)
            failed += not ex.pending
            pending += ex.pending
    if args.write:
        new = rewrite(lines, examples, args.pending)
        if new != lines:
            with open(args.page, 'w', encoding='utf-8') as f:
                f.write('\n'.join(new))
    print('%d examples: %d failing, %d pending examples differing from their intended '
          'output' % (len(examples), failed, pending))
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
