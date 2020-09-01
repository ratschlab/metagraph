#!/usr/bin/env python3

import argparse
import buckets
import cgi
import errno
import http.server
import json
import logging
import os
import socket
import subprocess
import time
import urllib
import urllib.parse

# TODO:
#  - maybe handle retries for pending jobs

args = None  # the parsed command line arguments

sra_generator = None  # class responsible for generating SRA ids to process
sra_download_gen = None  # the generator function in sra_generator that returns the next SRA id to be downloaded

# jobs done
downloaded_sras = set()  # the downloaded sras
built_sras = set()  # the sras for which a BOSS graph was generated
cleaned_sras = set()  # the sras for which the BOSS graph was cleaned (so processing is completed)
transferred_sras = set()  # the sras that were successfully transferred to permanent storage

# jobs currently being processed by a worker
pending_jobs = {}  # includes jobs that were pre-empted
pending_downloads = set()
pending_builds = set()
pending_cleans = set()
pending_transfers = set()

preempted_ids = set()

# set to true when all download ids were read
download_done = False

start_time = time.time()
op_time_ms = {}  # time spent for each operation, e.g. download, build, clean, transfer
total_download_size_MB = 0  # total MB of downloaded (and decompressed) data
total_build_size_MB = 0  # total MB of built data (excludes downloaded data that wasn't yet built)
total_clean_size_MB = 0  # total MB of cleaned data
total_transfer_size_MB = 0  # total MB of transferred data
sra_to_size = {}  # maps from an SRA id to its size

status_str = f"""
<html>
<head>
<title>Status of metagraph server</title>
</head>
<body>
<h3> Progress </h3>
<p>Processed: %s/%s</p>
<p>Uptime: %s</p>
<p>MB downloaded: %s</p>
<p>Download Time: %s (%s MB/s/machine)</p>
<p>MB built: %s</p>
<p>Build Time: %s (%s MB/s/machine)</p>
<p>MB Cleaned: %s</p>
<p>Clean Time: %s (%s MB/s/clean)</p>
<p>MB transferred: %s</p>
<p>Transfer Time: %s (%s MB/s/machine)</p>
<p>Overall processing speed (for all instances): %s MB/s</p>
<h3> Pending jobs </h3>
<p>Pending downloads: %s</p>
<p>Pending build: %s</p>
<p>Pending clean: %s</p>
<p>Pending transfer: %s</p>

<h3> Jobs preempted </h3>
<p>%s</p>

<h3> Jobs completed </h3>
<p>Completed downloads: %s (%s) </p>
<p>Completed build: %s (%s)</p>
<p>Completed clean: %s (%s)</p>
<p>Completed transfer: %s (%s)</p>

<h3> Download status </h3>
Download done: %s

</body>
</html>
"""


def internal_ip():
    try:
        return socket.gethostbyname(socket.gethostname())
    except socket.gaierror:
        return '127.0.0.1'  # this usually happens on dev laptops; cloud machines work fine


def publish_ip():
    with open('/tmp/server', 'w') as fp:
        fp.write(f'{internal_ip()}:{args.port}')
    if subprocess.call(['gsutil', 'cp', '/tmp/server', args.server_info], stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE) != 0:
        logging.error("Cannot publish server ip/port on Google Cloud Storage. Sorry, I tried.")
        exit(1)


def get_var(post_vars, var):
    var_list = post_vars.get(var.encode('ascii'))
    if not var_list:
        return None
    return var_list[0].decode('utf-8')


def millis_to_human(time_ms):
    seconds = int(time_ms) % 60
    minutes = int((time_ms / 60)) % 60
    hours = int((time_ms / (60 * 60))) % 24
    days = int((time_ms / (60 * 60 * 24)))
    return f'{days}d {hours}:{minutes}:{seconds}'


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
        global start_time
        uptime = millis_to_human(time.time() - start_time)
        download_time = millis_to_human(op_time_ms.get('download', 0))
        download_MBps = round(total_download_size_MB / op_time_ms.get('download', 1), 2)
        build_time = millis_to_human(op_time_ms.get('build', 0))
        build_MBps = round(total_build_size_MB / op_time_ms.get('build', 1), 2)
        clean_MBps = round(total_clean_size_MB / op_time_ms.get('clean', 1), 2)
        transfer_MBps = round(total_transfer_size_MB / op_time_ms.get('transfer', 1), 2)
        clean_time = millis_to_human(op_time_ms.get('clean', 0))
        transfer_time = millis_to_human(op_time_ms.get('transfer', 0))
        total_MBps = round(total_transfer_size_MB / (time.time() - start_time), 2)
        global sra_generator
        self.send_reply(200, status_str % (
            sra_generator.processed_sras, sra_generator.total_sras, uptime, round(total_download_size_MB, 2),
            download_time, download_MBps, total_build_size_MB, build_time, build_MBps, round(total_clean_size_MB, 2),
            clean_time, clean_MBps, round(total_transfer_size_MB, 2), transfer_time, transfer_MBps, total_MBps,
            pending_downloads, pending_builds, pending_cleans, pending_transfers, preempted_ids, len(downloaded_sras),
            downloaded_sras, len(built_sras), built_sras, len(cleaned_sras), cleaned_sras, len(transferred_sras),
            transferred_sras, download_done),
                        {'Content-type': 'text/html'})

    def handle_ack(self, operation, post_vars, add_sets, pending_operations):
        sra_id = get_var(post_vars, 'id')
        if None in (sra_id,):
            self.send_reply(400, 'id not specified')
            return False
        op_time_ms[operation] = op_time_ms.get(operation, 0) + int(get_var(post_vars, 'time'))
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
        if self.handle_ack('download', post_vars, [downloaded_sras, pending_builds], pending_downloads):
            global total_download_size_MB, sra_to_size
            size_MB = float(get_var(post_vars, 'size_mb'))
            if size_MB == 0:
                size_MB = float(get_var(post_vars, 'download_size_mb'))
            total_download_size_MB += size_MB
            sra_to_size[get_var(post_vars, 'id')] = size_MB

    def handle_ack_build(self, post_vars):
        if self.handle_ack('build', post_vars, [built_sras, pending_cleans], pending_builds):
            global total_build_size_MB, sra_to_size
            sra_id = get_var(post_vars, 'id')
            if sra_id in sra_to_size:
                total_build_size_MB += sra_to_size[sra_id]
            else:
                logging.warning(f'Sra {sra_id} was built but it\'s size can\'t be found. Server restart?')

    def handle_ack_clean(self, post_vars):
        if self.handle_ack('clean', post_vars, [cleaned_sras, pending_transfers], pending_cleans):
            # a transfer is automatically started after a clean operation by the client
            pending_transfers.add(get_var(post_vars, 'id'))
            global total_clean_size_MB, sra_to_size
            sra_id = get_var(post_vars, 'id')
            if sra_id in sra_to_size:
                total_clean_size_MB += sra_to_size[sra_id]
            else:
                logging.warning(f'Sra {sra_id} was cleaned but it\'s size can\'t be found. Server restart?')

    def handle_ack_transfer(self, post_vars):
        self.handle_ack('transfer', post_vars, [transferred_sras], pending_transfers)
        sra_id = get_var(post_vars, 'id')
        if sra_id in pending_jobs:
            del pending_jobs[sra_id]
            global total_transfer_size_MB, sra_to_size
            if sra_id in sra_to_size:
                total_transfer_size_MB += sra_to_size[sra_id]
            else:
                logging.warning(f'Sra {sra_id} was transferred but it\'s size can\'t be found. Server restart?')
        else:
            logging.warning(f'Cannot ack transfer of {sra_id}. Job was not pending. Server restart?')

    def handle_nack_download(self, post_vars):
        self.handle_nack('download', post_vars, pending_downloads)

    def handle_nack_build(self, post_vars):
        self.handle_nack('build', post_vars, pending_builds)

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
            if sra_id in pending_jobs:  # if server was restarted, pending_jobs info was lost :(
                preempted_ids.add(sra_id)
            pending_downloads.discard(sra_id)
            pending_builds.discard(sra_id)
            pending_cleans.discard(sra_id)
            pending_transfers.discard(sra_id)
            downloaded_sras.discard(sra_id)
            built_sras.discard(sra_id)
            cleaned_sras.discard(sra_id)
        self.send_reply(200, "")

    def send_reply(self, code, message, headers={}):
        self.send_response(code)
        for k, v in headers.items():
            self.send_header(k, v)
        self.end_headers()
        self.flush_headers()
        self.wfile.write(message.encode('utf-8'))
        self.wfile.flush()

    def get_post_vars(self, query_string):
        ctype, pdict = cgi.parse_header(self.headers['content-type'])
        if ctype != 'application/x-www-form-urlencoded':
            self.send_reply(400, "Bad content-type, only application/x-www-form-urlencoded accepted")
            return None
        return urllib.parse.parse_qs(query_string, keep_blank_values=True)

    def do_GET(self):
        logging.info(f'GET {self.path}')
        parsed_url = urllib.parse.urlparse(self.path)
        if parsed_url.path == '/jobs':
            self.handle_get_job()
        elif parsed_url.path == '/status':
            self.handle_get_status()
        else:
            self.send_reply(404, f'Invalid path: {self.path}\n')

    def do_POST(self):
        parsed_url = urllib.parse.urlparse(self.path)
        length = int(self.headers.get('content-length', 0))
        query_string = self.rfile.read(length)
        post_vars = self.get_post_vars(query_string)
        logging.info(f'POST {parsed_url.path} {query_string.decode("utf-8")}')
        if not post_vars:
            self.send_reply(400, f'No POST parameters specified.')
            return
        if parsed_url.path == '/jobs/ack/download':
            self.handle_ack_download(post_vars)
        elif parsed_url.path == '/jobs/ack/build':
            self.handle_ack_build(post_vars)
        elif parsed_url.path == '/jobs/ack/clean':
            self.handle_ack_clean(post_vars)
        elif parsed_url.path == '/jobs/ack/transfer':
            self.handle_ack_transfer(post_vars)
        elif parsed_url.path == '/jobs/nack/download':
            self.handle_nack_download(post_vars)
        elif parsed_url.path == '/jobs/nack/build':
            self.handle_nack_build(post_vars)
        elif parsed_url.path == '/jobs/nack/clean':
            self.handle_nack_clean(post_vars)
        elif parsed_url.path == '/jobs/nack/transfer':
            self.handle_nack_transfer(post_vars)
        elif parsed_url.path == '/jobs/preempt':
            self.handle_preempt(post_vars)
        else:
            self.send_reply(404, f'Invalid path: {self.path}\n')


def load_file_set(filename):
    result = set()
    try:
        with open(filename) as fp:
            for line in fp:
                if len(line.rstrip()) == 0:
                    continue
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
    global downloaded_sras, built_sras, cleaned_sras, transferred_sras
    state_map = {'downloaded': load_file_set(filename),
                 'built': load_file_set(os.path.join(args.output_dir, 'succeed_build.id')),
                 'cleaned': load_file_set(os.path.join(args.output_dir, 'succeed_clean.id')),
                 'transferred': load_file_set(os.path.join(args.output_dir, 'succeed_transfer.id')),
                 '!downloaded': load_file_set(os.path.join(args.output_dir, 'failed_download.id')),
                 '!built': load_file_set(os.path.join(args.output_dir, 'failed_build.id')),
                 '!cleaned': load_file_set(os.path.join(args.output_dir, 'failed_clean.id')),
                 '!transferred': load_file_set(os.path.join(args.output_dir, 'failed_transfer.id'))}

    checkpoint = args.checkpoint.split('|')
    to_eliminate = set()
    for state in checkpoint:
        state_size = len(state_map[state])
        if state_size != 0:
            logging.info(f'Eliminating {state} of size {state_size}')
        to_eliminate.update(state_map[state])
    global sra_download_gen, sra_generator
    sra_generator = buckets.Sra(data_files, args.add_gcloud_bucket, args.data_dir, to_eliminate)
    sra_download_gen = sra_generator.next_item()


def init_logging():
    args.output_dir = os.path.expanduser(args.output_dir)
    if not os.path.exists(args.output_dir):
        try:
            os.makedirs(args.output_dir)
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    logging.basicConfig(format='%(asctime)s %(message)s', level=logging.DEBUG, datefmt='%Y-%m-%d %H:%M:%S')
    file_handler = logging.FileHandler(f'{args.output_dir}/server.log')
    file_handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s %(message)s')
    file_handler.setFormatter(formatter)
    logging.getLogger().addHandler(file_handler)


def check_env():
    """ Make sure all the necessary software is in place to successfully run the serer """
    if subprocess.call(['./prereq.sh']) != 0:
        logging.error("Some prerequisites are missing on this machine. Bailing out.")
        exit(1)


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--port', default=8000, type=int, help='HTTP Port on which the server runs')
    parser.add_argument(
        '--data_dir',
        default=os.path.expanduser('~/Downloads/sra_test/'),
        help='Location of the directory containing the input data')
    parser.add_argument(
        '--output_dir',
        default=os.path.expanduser('~/.metagraph/'),
        help='Location of the directory containing the input data')
    parser.add_argument('--add_gcloud_bucket', default=True, dest='add_gcloud_bucket', action='store_true',
                        help='Whether to add the NCBI gcloud bucket id for each sra')
    parser.add_argument('--noadd_gcloud_bucket', default=True, dest='add_gcloud_bucket', action='store_false',
                        help='Whether to NOT add the NCBI gcloud bucket id for each sra')
    parser.set_defaults(feature=True)
    # leave this to 1 to avoid race conditions, or fix race conditions :)
    parser.add_argument('--worker_count', default=1, help='Number of workers processing data')
    parser.add_argument('--checkpoint', default='transferred|!downloaded|!built|!cleaned',
                        help='Which checkpointed SRAs to eliminate from processing, '
                             'a combination of downloaded/built/cleaned/transferred optionally preceded by !')
    parser.add_argument('--server_info', default=None,
                        help='Where to publish the server host/port on gcs')

    global args
    args = parser.parse_args()


if __name__ == '__main__':
    parse_args()
    init_logging()
    check_env()
    init_state()
    logging.info(f'Starting server on port {args.port}...')
    httpd = http.server.HTTPServer(('', args.port), SimpleHTTPRequestHandler)
    logging.info(f'Publishing server address {internal_ip()}:{args.port}...')
    publish_ip()

    httpd.serve_forever()
