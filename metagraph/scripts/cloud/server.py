import argparse
import cgi
import errno
import http.server
import json
import os
import urllib

args = None  # the parsed command line arguments

downloaded_sras = map()  # the downloaded sras
created_sras = map()  # the sras for which a BOSS graph was generated
cleaned_sras = map()  # the sras for which the BOSS graph was cleaned (so processing is completed)
transferred_sras = map() # the sras that were successfully transferred to permanent storage

sra_download_gen = None  # the generator function that returns the next SRA id to be downloaded
to_create_sras = map()
to_clean_sras = map()
to_transfer_sras = map()

pending_downloads = set()
pending_creates = set()
pending_cleans = set()
pending_transfers = set()

class Sra:
    def __init__(self, data_files):
        self.data_files = data_files

    def next_item(self):
        """Lazy function (generator) to read a sequence of files line by line."""
        for file in data_files:
            with open(file) as fp:
                for line in fp:
                    if line in downloaded_sras:  # already processed
                        continue
                    yield line
        return None

class SimpleHTTPRequestHandler(http.server.BaseHTTPRequestHandler):
    def handle_get_job(self, params):
        try:
            download_jobs = params.get('download')
            create_jobs = params.get('create')
            clean_jobs = params.get('clean')
            if None in (download_jobs, create_jobs, clean_jobs):
                self.send_reply(400, "One of 'download', 'create' or 'clean' parameters is missing\n")
                return
            response = {}
            print(f'{download_jobs} {create_jobs}')
            if download_jobs[0] == '0':
                sra_id = next(sra_download_gen)
                response['download'] = sra_id
                pending_downloads.add(sra_id)
            create_jobs_int = int(create_jobs)
            clean_jobs_int = int(clean_jobs)
            if clean_jobs_int > 0 and clean_jobs_int < 4:


            self.send_header("Content-type", "application/json")
            self.send_reply(200, json.dumps(response))
        except StopIteration:
            self.send_reply(204, "All done!")  # no content
        self.end_headers()

    def handle_ack_download(self, postvars):
        sra_id = postvars.get(b'sra_id')
        location = postvars.get(b'location')
        if None in (sra_id, location):
            self.send_reply(400, 'sra_id or location not specified')
            return
        filename = os.path.join(args.output_dir, 'downloaded_sras')
        sra_id_str = sra_id[0].decode('utf-8')
        location_str = location[0].decode('utf-8')
        with open(filename, 'a') as fp:
            fp.write(f'{sra_id_str} {location_str}\n')
        downloaded_sras.add({sra_id_str: location_str})
        self.send_reply(200, "")

        pending_downloads.remove(sra_id_str)

    def handle_ack_create(self, postvars):
        sra_id = postvars.get(b'sra_id')
        location = postvars.get(b'location')
        if None in (sra_id, location):
            self.send_reply(400, 'sra_id or location not specified')
            return
        filename = os.path.join(args.output_dir, 'created_sras')
        sra_id_str = sra_id[0].decode('utf-8')
        location_str = location[0].decode('utf-8')
        with open(filename, 'a') as fp:
            fp.write(f'{sra_id_str} {location_str}\n')
        downloaded_sras.add({sra_id_str: location_str})
        self.send_reply(200, "")

        pending_creates.remove(sra_id_str)

    def do_GET(self):
        parsed_url = urllib.parse.urlparse(self.path)
        params = urllib.parse.parse_qs(parsed_url.query)
        if parsed_url.path == '/jobs':
            self.handle_get_job(params)
        else:
            self.send_reply(404, f'Invalid path: {self.path}')

    def send_reply(self, code, message):
        self.send_response(code)
        self.end_headers()
        self.wfile.write(message.encode('utf-8'))

    def get_postvars(self):
        ctype, pdict = cgi.parse_header(self.headers['content-type'])
        if ctype != 'application/x-www-form-urlencoded':
            self.send_reply(400, "Bad content-type, only x-www-form-urlencoded accepted")
            return None
        length = int(self.headers['content-length'])
        return urllib.parse_qs(self.rfile.read(length), keep_blank_values=True)

    def do_POST(self):
        parsed_url = urllib.parse.urlparse(self.path)
        postvars = self.get_postvars()
        if not postvars:
            return
        if parsed_url.path == '/jobs/ack/download':
            self.handle_ack_download(postvars)
        if parsed_url.path == '/jobs/ack/create':
            self.handle_ack_create(postvars)
        if parsed_url.path == '/jobs/ack/clean':
            self.handle_ack_clean(postvars)
        if parsed_url.path == '/jobs/ack/transfer':
            self.handle_ack_transfer(postvars)
        else:
            self.send_reply(404, f'Invalid path: {self.path}')


def init_state():
    filename = os.path.join(args.output_dir, 'downloaded_sras')
    if not os.path.exists(os.path.dirname(filename)):
        try:
            os.makedirs(os.path.dirname(filename))
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    try:
        with open(filename) as fp:
            for line in fp:
                downloaded_sras.add(line)
    except FileNotFoundError:  # no downloaded sras yet, that's fine
        return


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

    sra_download_gen = Sra(data_files).next_item()
    sra_create_gen = Sra(data_files).next_item()
    sra_clean_gen = Sra(data_files).next_item()

    httpd = http.server.HTTPServer(('localhost', args.port), SimpleHTTPRequestHandler)
    print(f'Starting server on port {args.port}')
    httpd.serve_forever()
