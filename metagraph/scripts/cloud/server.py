from http.server import HTTPServer, BaseHTTPRequestHandler
import argparse
import errno
import os
from cgi import parse_header
from urllib.parse import parse_qs

acked_sras = set()  # the acked sras
args = None  # the parsed command line arguments
sra_gen = None  # the generator function that returns the next SRA id to be processed

class Sra:
    def __init__(self, data_files):
        self.data_files = data_files

    def next_item(self):
        """Lazy function (generator) to read a sequence of files line by line."""
        for file in data_files:
            with open(file) as fp:
                for line in fp:
                    if line in acked_sras:  # already processed
                        continue
                    yield line
        return None


class SimpleHTTPRequestHandler(BaseHTTPRequestHandler):
    def do_GET(self):
        if self.path == '/jobs':
            try:
                self.send_header("Content-type", "application/json")
                self.send_reply(200, next(sra_gen))
            except StopIteration:
                self.send_reply(204, "All done!")  # no content
            self.end_headers()
        else:
            self.send_reply(404, f'Invalid path: {self.path}')

    def send_reply(self, code, message):
        self.send_response(code)
        self.end_headers()
        self.wfile.write(message.encode('utf-8'))

    def get_postvars(self):
        ctype, pdict = parse_header(self.headers['content-type'])
        if ctype != 'application/x-www-form-urlencoded':
            self.send_reply(400, "Bad content-type, only x-www-form-urlencoded accepted")
            return None
        length = int(self.headers['content-length'])
        return parse_qs(self.rfile.read(length), keep_blank_values=1)

    def do_POST(self):
        if self.path == '/jobs/ack':
            postvars = self.get_postvars()
            if not postvars:
                return
            sra_id = postvars.get(b'sra_id')
            if not sra_id:
                self.send_reply(400, 'sra_id is not specified')
                return
            filename = os.path.join(args.output_dir, 'acked_sras')
            if not os.path.exists(os.path.dirname(filename)):
                try:
                    os.makedirs(os.path.dirname(filename))
                except OSError as exc:  # Guard against race condition
                    if exc.errno != errno.EEXIST:
                        raise
            with open(filename, 'a') as fp:
                fp.write(f'{sra_id[0].decode("utf-8")}\n')
            acked_sras.add(sra_id[0])
            self.send_reply(200, "")
        else:
            self.send_reply(404, f'Invalid path: {self.path}')


def init_state():
    filename = os.path.join(args.output_dir, 'acked_sras')
    with open(filename) as fp:
        for line in fp:
            acked_sras.add(line)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--port', default=8000, help='HTTP Port on which the server runs')
    parser.add_argument(
        '--data_dir',
        default=os.path.expanduser('~/Downloads/sra_test/'),
        help='Location of the directory containing the input data')
    parser.add_argument(
        '--output_dir',
        default=os.path.expanduser('~/.metagraph/'),
        help='Location of the directory containing the input data')

    args = parser.parse_args()

    data_files = [os.path.join(args.data_dir, f) for f in os.listdir(args.data_dir) if
                  os.path.isfile(os.path.join(args.data_dir, f))]

    if len(data_files) == 0:
        print(f"No files found in '{args.data_dir}'. Exiting.")
        exit(1)

    init_state()

    sra_gen = Sra(data_files).next_item()

    httpd = HTTPServer(('localhost', args.port), SimpleHTTPRequestHandler)
    print(f'Starting server on port {args.port}')
    httpd.serve_forever()
