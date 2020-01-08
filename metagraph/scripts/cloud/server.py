import argparse
import cgi
import errno
import http.server
import json
import logging
import os
import urllib

# TODO:
# - handle retries for pending jobs

args = None  # the parsed command line arguments

# jobs done
downloaded_sras = {}  # the downloaded sras
created_sras = {}  # the sras for which a BOSS graph was generated
cleaned_sras = {}  # the sras for which the BOSS graph was cleaned (so processing is completed)
transferred_sras = {}  # the sras that were successfully transferred to permanent storage

# jobs in the queue (waiting for a worker)
sra_download_gen = None  # the generator function that returns the next SRA id to be downloaded
to_create_sras = {}
to_clean_sras = {}
to_transfer_sras = {}

# jobs currently being processed by a worker
pending_downloads = set()
pending_creates = set()
pending_cleans = set()
pending_transfers = set()

# set to true when all download ids were read
download_done = False

logger = logging.getLogger('metagraph-server')

status_str = f"""
<html>
<head>
<title>Status of metagraph server</title>
</head>
<body>
<h3> Pending jobs </h3>
<p>Pending downloads: %s</p>
<p>Pending create: %s</p>
<p>Pending clean: %s</p>
<p>Pending transfer: %s</p>

<h3> Jobs waiting to be scheduled </h3>
<p>Waiting create: %s</p>
<p>Waiting clean: %s</p>
<p>Waiting transfer: %s</p>

<h3> Jobs completed </h3>
<p>Completed downloads: %s </p>
<p>Completed create: %s</p>
<p>Completed clean: %s</p>
<p>Completed transfer: %s</p>

</body>
</html>
"""


class Sra:
    """ Reads a list of files containing SRA ids and returns the ids one by one for further processing """

    def __init__(self, data_files):
        self.data_files = data_files

    def next_item(self):
        """Lazy function (generator) to read a sequence of files line by line."""
        for file in data_files:
            with open(file) as fp:
                for line in fp:
                    line = line.rstrip()
                    if line in downloaded_sras:  # already processed
                        print(f'Lone {line} already downloaded')
                        continue
                    yield line
        return None


class SimpleHTTPRequestHandler(http.server.BaseHTTPRequestHandler):
    """ Processes requests for new jobs or for acknowledging finished jobs """

    def handle_get_job(self, params):
        """ Handles a request for a new job """

        download_jobs = params.get('download')
        create_jobs = params.get('create')
        clean_jobs = params.get('clean')
        if None in (download_jobs, create_jobs, clean_jobs):
            self.send_reply(400, "One of 'download', 'create' or 'clean' parameters is missing\n")
            return
        global download_done
        if download_done and not (
                to_create_sras or to_clean_sras or to_transfer_sras or pending_downloads or pending_creates
                or pending_cleans or pending_transfers):
            self.send_reply(204, 'All done, feel free to exit and pat yourself on the back')
        response = {}
        if download_jobs[0] == '0':
            try:
                sra_id = next(sra_download_gen)
                response['download'] = {'id': sra_id}
                pending_downloads.add(sra_id)
            except StopIteration:
                download_done = True  # nothing else to download
        create_jobs_int = int(create_jobs[0])
        clean_jobs_int = int(clean_jobs[0])
        if not to_create_sras or (0 < clean_jobs_int < 4):  # there is room for more clean jobs
            if to_clean_sras:
                sra_id, dbg_file = to_clean_sras.popitem()
                response['clean'] = {'id': sra_id, 'location': dbg_file}
                pending_cleans.add(sra_id)
        elif create_jobs_int == 0:
            if to_create_sras:
                sra_id, directory = to_create_sras.popitem()
                response['create'] = {'id': sra_id, 'location': directory}
                pending_creates.add(sra_id)
        self.send_reply(200, json.dumps(response), {'Content-type': 'application/json'})

    def handle_get_status(self):
        self.send_reply(200, status_str % (
            pending_downloads, pending_creates, pending_cleans, pending_transfers, to_create_sras, to_clean_sras,
            to_transfer_sras, downloaded_sras, created_sras, cleaned_sras, transferred_sras),
                        {'Content-type': 'text/html'})

    def handle_ack(self, operation, post_vars, add_maps, remove_map):
        sra_id = post_vars.get(b'id')
        location = post_vars.get(b'location')
        if None in (sra_id, location):
            self.send_reply(400, 'id or location not specified')
            return
        filename = os.path.join(args.output_dir, f'{operation}ed_sras')
        sra_id_str = sra_id[0].decode('utf-8')
        location_str = location[0].decode('utf-8')
        with open(filename, 'a') as fp:
            fp.write(f'{sra_id_str} {location_str}\n')
        for m in add_maps:
            m[sra_id_str] = location_str
        self.send_reply(200, "")
        try:
            remove_map.remove(sra_id_str)
        except KeyError:
            logger.warning(f'Acknowledging nonexistent pending {operation} {sra_id_str}')

    def handle_nack(self, operation, post_vars, remove_map):
        sra_id = post_vars.get(b'id')
        if sra_id is None:
            self.send_reply(400, 'id not specified')
            return
        filename = os.path.join(args.output_dir, f'failed_{operation}.id')
        sra_id_str = sra_id[0].decode('utf-8')
        with open(filename, 'a') as fp:
            fp.write(f'{sra_id_str}\n')
        self.send_reply(200, "")
        try:
            remove_map.remove(sra_id_str)
        except KeyError:
            logger.warning(f'Acknowledging nonexistent pending {operation} {sra_id_str}')

    def handle_ack_download(self, postvars):
        self.handle_ack('download', postvars, [downloaded_sras, to_create_sras], pending_downloads)

    def handle_ack_create(self, postvars):
        self.handle_ack('create', postvars, [created_sras, to_clean_sras], pending_creates)

    def handle_ack_clean(self, postvars):
        self.handle_ack('clean', postvars, [cleaned_sras, to_transfer_sras], pending_cleans)

    def handle_nack_download(self, postvars):
        self.handle_nack('download', postvars, pending_downloads)

    def handle_nack_create(self, postvars):
        self.handle_nack('create', postvars, pending_creates)

    def handle_nack_clean(self, postvars):
        self.handle_nack('clean', postvars, pending_cleans)

    def send_reply(self, code, message, headers={}):
        self.send_response(code)
        for k, v in headers.items():
            self.send_header(k, v)
        self.end_headers()
        self.flush_headers()
        self.wfile.write(message.encode('utf-8'))
        self.wfile.flush()

    def get_postvars(self):
        ctype, pdict = cgi.parse_header(self.headers['content-type'])
        if ctype != 'application/x-www-form-urlencoded':
            self.send_reply(400, "Bad content-type, only x-www-form-urlencoded accepted")
            return None
        length = int(self.headers['content-length'])
        return urllib.parse.parse_qs(self.rfile.read(length), keep_blank_values=True)

    def do_GET(self):
        parsed_url = urllib.parse.urlparse(self.path)
        params = urllib.parse.parse_qs(parsed_url.query)
        if parsed_url.path == '/jobs':
            self.handle_get_job(params)
        elif parsed_url.path == '/status':
            self.handle_get_status()
        else:
            self.send_reply(404, f'Invalid path: {self.path}\n')

    def do_POST(self):
        parsed_url = urllib.parse.urlparse(self.path)
        postvars = self.get_postvars()
        if not postvars:
            return
        if parsed_url.path == '/jobs/ack/download':
            self.handle_ack_download(postvars)
        elif parsed_url.path == '/jobs/ack/create':
            self.handle_ack_create(postvars)
        elif parsed_url.path == '/jobs/ack/clean':
            self.handle_ack_clean(postvars)
        if parsed_url.path == '/jobs/nack/download':
            self.handle_nack_download(postvars)
        elif parsed_url.path == '/jobs/nack/create':
            self.handle_nack_create(postvars)
        elif parsed_url.path == '/jobs/nack/clean':
            self.handle_nack_clean(postvars)
        else:
            self.send_reply(404, f'Invalid path: {self.path}\n')


def load_file_dict(filename):
    result = {}
    try:
        with open(filename) as fp:
            for line in fp:
                tokens = line.split()
                assert len(tokens) == 2, f'Invalid line in downloads sras {line}'
                result[tokens[0]] = tokens[1]
    except FileNotFoundError:  # no downloaded sras yet, that's fine
        return result
    return result


def init_state():
    filename = os.path.join(args.output_dir, 'downloaded_sras')
    if not os.path.exists(os.path.dirname(filename)):
        try:
            os.makedirs(os.path.dirname(filename))
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    global downloaded_sras, created_sras, cleaned_sras, transferred_sras, to_create_sras, to_clean_sras, to_transfer_sras
    downloaded_sras = load_file_dict(filename)
    created_sras = load_file_dict(os.path.join(args.output_dir, 'created_sras'))
    cleaned_sras = load_file_dict(os.path.join(args.output_dir, 'cleaned_sras'))
    transferred_sras = load_file_dict(os.path.join(args.output_dir, 'transferred_sras'))

    to_create_sras = {k: downloaded_sras[k] for k in set(downloaded_sras) - set(created_sras)}
    to_clean_sras = {k: created_sras[k] for k in set(created_sras) - set(cleaned_sras)}
    to_transfer_sras = {k: cleaned_sras[k] for k in set(cleaned_sras) - set(transferred_sras)}


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
                  os.path.isfile(os.path.join(args.data_dir, f)) and f.endswith('ids')]

    if len(data_files) == 0:
        print(f"No files found in '{args.data_dir}'. Exiting.")
        exit(1)

    init_state()

    sra_download_gen = Sra(data_files).next_item()

    httpd = http.server.HTTPServer(('localhost', args.port), SimpleHTTPRequestHandler)
    print(f'Starting server on port {args.port}')
    httpd.serve_forever()
