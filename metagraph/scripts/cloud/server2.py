#!/usr/bin/env python3

import argparse
import buckets
import cgi
import errno
import http.server
import json
import logging
import os
import urllib

# TODO:
#  - maybe handle retries for pending jobs

args = None  # the parsed command line arguments

sra_download_gen = None  # the generator function that returns the next SRA id to be downloaded

# jobs done
downloaded_sras = set()  # the downloaded sras
created_sras = set()  # the sras for which a BOSS graph was generated
cleaned_sras = set()  # the sras for which the BOSS graph was cleaned (so processing is completed)
transferred_sras = set()  # the sras that were successfully transferred to permanent storage

# jobs currently being processed by a worker
pending_jobs = {}  # includes jobs that were pre-empted
pending_downloads = set()
pending_creates = set()
pending_cleans = set()
pending_transfers = set()

preempted_ids = set()

# set to true when all download ids were read
download_done = False

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

<h3> Jobs preempted </h3>
<p>%s</p>

<h3> Jobs completed </h3>
<p>Completed downloads: %s </p>
<p>Completed create: %s</p>
<p>Completed clean: %s</p>
<p>Completed transfer: %s</p>

<h3> Download status </h3>
Download done: %s

</body>
</html>
"""


def get_var(post_vars, var):
    var_list = post_vars.get(var.encode('ascii'))
    if not var_list:
        return None
    return var_list[0].decode('utf-8')


class SimpleHTTPRequestHandler(http.server.BaseHTTPRequestHandler):
    """ Processes requests for new jobs or for acknowledging finished jobs """

    def handle_get_job(self):
        """ Handles a request for a new job """

        global download_done
        response = {}
        try:
            if preempted_ids:
                sra_id = preempted_ids.pop()
                bucket = pending_jobs[sra_id]
            else:
                sra_id_line = next(sra_download_gen)
                sra_id = sra_id_line[0]
                if len(sra_id_line) == 2:
                    bucket = sra_id_line[1]
                else:
                    bucket = None
                pending_jobs[sra_id] = bucket
            response['download'] = {'id': sra_id}
            if bucket:
                response['download']['bucket'] = bucket
            pending_downloads.add(sra_id)

        except StopIteration:
            download_done = True  # nothing else to download

        self.send_reply(200, json.dumps(response), {'Content-type': 'application/json'})

    def handle_get_status(self):
        self.send_reply(200, status_str % (
            pending_downloads, pending_creates, pending_cleans, pending_transfers, preempted_ids, downloaded_sras,
            created_sras, cleaned_sras, transferred_sras, download_done),
                        {'Content-type': 'text/html'})

    def handle_ack(self, operation, post_vars, add_sets, pending_operations):
        sra_id = get_var(post_vars, 'id')
        location = get_var(post_vars, 'location')
        if None in (sra_id, location):
            self.send_reply(400, 'id or location not specified')
            return False
        filename = os.path.join(args.output_dir, f'succeed_{operation}.id')
        with open(filename, 'a') as fp:
            fp.write(f'{sra_id}\n')
        for s in add_sets:
            s.add(sra_id)
        self.send_reply(200, "")
        try:
            pending_operations.remove(sra_id)
        except KeyError:
            logging.warning(f'Acknowledging nonexistent pending {operation} {sra_id}')
        return True

    def handle_nack(self, operation, post_vars, pending_operations):
        sra_id = get_var(post_vars, 'id')
        if sra_id is None:
            self.send_reply(400, 'id not specified')
            return
        filename = os.path.join(args.output_dir, f'failed_{operation}.id')
        with open(filename, 'a') as fp:
            fp.write(f'{sra_id}\n')
        self.send_reply(200, "")
        try:
            pending_operations.remove(sra_id)
        except KeyError:
            logging.warning(f'Acknowledging nonexistent pending {operation} {sra_id}')

    def handle_ack_download(self, post_vars):
        self.handle_ack('download', post_vars, [downloaded_sras, pending_creates], pending_downloads)

    def handle_ack_create(self, post_vars):
        self.handle_ack('create', post_vars, [created_sras, pending_cleans], pending_creates)

    def handle_ack_clean(self, post_vars):
        if self.handle_ack('clean', post_vars, [cleaned_sras, pending_transfers], pending_cleans):
            # a transfer is automatically started after a clean operation by the client
            pending_transfers.add(get_var(post_vars, 'id'))

    def handle_ack_transfer(self, post_vars):
        self.handle_ack('transfer', post_vars, [transferred_sras], pending_transfers)
        sra_id = get_var(post_vars, 'id')
        del pending_jobs[sra_id]

    def handle_nack_download(self, post_vars):
        self.handle_nack('download', post_vars, pending_downloads)

    def handle_nack_create(self, post_vars):
        self.handle_nack('create', post_vars, pending_creates)

    def handle_nack_clean(self, post_vars):
        self.handle_nack('clean', post_vars, pending_cleans)

    def handle_nack_transfer(self, post_vars):
        self.handle_nack('transfer', post_vars, pending_transfers)

    def handle_preempt(self, post_vars):
        sra_ids = get_var(post_vars, 'ids')
        if sra_ids is None:
            self.send_reply(400, 'No ids specified. Need at least one')
            return
        for sra_id in sra_ids.split(','):
            sra_id = sra_id.strip()
            preempted_ids.add(sra_id)
            pending_downloads.remove(sra_id)
            pending_creates.remove(sra_id)
            pending_cleans.remove(sra_id)
            pending_transfers.remove(sra_id)

    def send_reply(self, code, message, headers={}):
        self.send_response(code)
        for k, v in headers.items():
            self.send_header(k, v)
        self.end_headers()
        self.flush_headers()
        self.wfile.write(message.encode('utf-8'))
        self.wfile.flush()

    def get_post_vars(self):
        ctype, pdict = cgi.parse_header(self.headers['content-type'])
        if ctype != 'application/x-www-form-urlencoded':
            self.send_reply(400, "Bad content-type, only application/x-www-form-urlencoded accepted")
            return None
        length = int(self.headers['content-length'])
        return urllib.parse.parse_qs(self.rfile.read(length), keep_blank_values=True)

    def do_GET(self):
        parsed_url = urllib.parse.urlparse(self.path)
        if parsed_url.path == '/jobs':
            self.handle_get_job()
        elif parsed_url.path == '/status':
            self.handle_get_status()
        else:
            self.send_reply(404, f'Invalid path: {self.path}\n')

    def do_POST(self):
        parsed_url = urllib.parse.urlparse(self.path)
        post_vars = self.get_post_vars()
        logging.debug(f'POST {parsed_url} {post_vars}')
        if not post_vars:
            self.send_reply(400, f'No POST parameters specified.')
            return
        if parsed_url.path == '/jobs/ack/download':
            self.handle_ack_download(post_vars)
        elif parsed_url.path == '/jobs/ack/create':
            self.handle_ack_create(post_vars)
        elif parsed_url.path == '/jobs/ack/clean':
            self.handle_ack_clean(post_vars)
        elif parsed_url.path == '/jobs/ack/transfer':
            self.handle_ack_transfer(post_vars)
        elif parsed_url.path == '/jobs/nack/download':
            self.handle_nack_download(post_vars)
        elif parsed_url.path == '/jobs/nack/create':
            self.handle_nack_create(post_vars)
        elif parsed_url.path == '/jobs/nack/clean':
            self.handle_nack_clean(post_vars)
        elif parsed_url.path == '/jobs/preempt':
            self.handle_preempt(post_vars)
        else:
            self.send_reply(404, f'Invalid path: {self.path}\n')


def load_file_set(filename):
    result = set()
    try:
        with open(filename) as fp:
            for line in fp:
                tokens = line.split()
                assert len(tokens) == 1, f'Invalid line in downloads sras {line}'
                result.add(tokens[0])
    except FileNotFoundError:  # no downloaded sras yet, that's fine
        return result
    return result


def init_state():
    args.data_dir = os.path.expanduser(args.data_dir)
    if os.path.isfile(args.data_dir):
        logging.error(
            f' --data_dir must point to a directory (containing files with SRA ids), not to a file (e.g. {os.path.dirname(args.data_dir)} instead of {args.data_dir})')
        exit(1)
    if not os.path.exists(args.data_dir):
        logging.error(
            f' Directory {args.data_dir} does not exists. Please specify a valid input directory via the --data_dir flag')
        exit(1)
    data_files = [os.path.join(args.data_dir, f) for f in os.listdir(args.data_dir) if
                  os.path.isfile(os.path.join(args.data_dir, f)) and f.endswith('ids')]

    if len(data_files) == 0:
        logging.fatal(f"No files found in '{args.data_dir}'. Exiting.")
        exit(1)

    filename = os.path.join(args.output_dir, 'succeed_download.id')
    if not os.path.exists(os.path.dirname(filename)):
        try:
            os.makedirs(os.path.dirname(filename))
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    global downloaded_sras, created_sras, cleaned_sras, transferred_sras
    downloaded_sras = load_file_set(filename)
    created_sras = load_file_set(os.path.join(args.output_dir, 'succeed_create.id'))
    cleaned_sras = load_file_set(os.path.join(args.output_dir, 'succeed_clean.id'))
    transferred_sras = load_file_set(os.path.join(args.output_dir, 'succeed_transfer.id'))

    global sra_download_gen
    sra_download_gen = buckets.Sra(data_files, args.add_gcloud_bucket, args.data_dir, downloaded_sras).next_item()


def init_logging():
    args.output_dir = os.path.expanduser(args.output_dir)
    if not os.path.exists(args.output_dir):
        try:
            os.makedirs(args.output_dir)
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    logging.basicConfig(format='%(asctime)s %(message)s', level=logging.DEBUG)
    file_handler = logging.FileHandler("{0}/{1}.log".format(args.output_dir, 'server'))
    file_handler.setLevel(logging.DEBUG)
    logging.getLogger().addHandler(file_handler)


def parse_args():
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
    parser.add_argument('--add_gcloud_bucket', default=True,
                        help='Whether to add the NCBI gcloud bucket id for each sra')
    parser.add_argument('--worker_count', default=1, help='Number of workers processing data')

    global args
    args = parser.parse_args()


if __name__ == '__main__':
    parse_args()

    init_logging()

    init_state()

    httpd = http.server.HTTPServer(('', args.port), SimpleHTTPRequestHandler)
    logging.info(f'Starting server on port {args.port}')
    httpd.serve_forever()
