#!/usr/bin/env python3

import argparse
import buckets
import cgi
import errno
import fileinput
import http.server
import json
import logging
import os
import urllib

# TODO:
#  - maybe handle retries for pending jobs

args = None  # the parsed command line arguments

# jobs done
downloaded_sras = {}  # the downloaded sras
built_sras = {}  # the sras for which a BOSS graph was generated
cleaned_sras = {}  # the sras for which the BOSS graph was cleaned (so processing is completed)
transferred_sras = {}  # the sras that were successfully transferred to permanent storage

# jobs in the queue (waiting for a worker)
sra_download_gen = None  # the generator function that returns the next SRA id to be downloaded
to_build_sras = {}
to_clean_sras = {}
to_transfer_sras = {}

# jobs currently being processed by a worker
pending_downloads = set()
pending_builds = set()
pending_cleans = set()
pending_transfers = set()

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
<p>Pending build: %s</p>
<p>Pending clean: %s</p>
<p>Pending transfer: %s</p>

<h3> Jobs waiting to be scheduled </h3>
<p>Waiting build: %s</p>
<p>Waiting clean: %s</p>
<p>Waiting transfer: %s</p>

<h3> Jobs completed </h3>
<p>Completed downloads: %s </p>
<p>Completed build: %s</p>
<p>Completed clean: %s</p>
<p>Completed transfer: %s</p>

</body>
</html>
"""


def get_var(post_vars, var):
    var_list = post_vars.get(var.encode('ascii'))
    if not var_list:
        return None
    return var_list[0].decode('utf-8')


def should_download():
    """ Returns true if it makes sense to download more data for processing, e.g. if we don't already have a
    large number of downloads waiting to be processed"""

    return len(to_build_sras) < 10 and len(to_build_sras) < 3 * args.worker_count


def should_build():
    """ Returns true if it makes sense to build rather than clean """

    return len(to_build_sras) <= len(to_clean_sras) or len(to_clean_sras) < 4 * args.worker_count


class SimpleHTTPRequestHandler(http.server.BaseHTTPRequestHandler):
    """ Processes requests for new jobs or for acknowledging finished jobs """

    def handle_get_job(self, params):
        """ Handles a request for a new job """

        download_jobs = params.get('download')
        build_jobs = params.get('build')
        clean_jobs = params.get('clean')
        if None in (download_jobs, build_jobs, clean_jobs):
            self.send_reply(400, "One of 'download', 'build' or 'clean' parameters is missing\n")
            return
        global download_done
        if download_done and not (
                to_build_sras or to_clean_sras or to_transfer_sras or pending_downloads or pending_builds
                or pending_cleans or pending_transfers):
            self.send_reply(204, 'All done, feel free to exit and pat yourself on the back')
        response = {}
        if download_jobs[0] == '0' and should_download():
            try:
                sra_id_line = next(sra_download_gen)
                response['download'] = {'id': sra_id_line[0]}
                if len(sra_id_line) == 2:
                    response['download']['bucket'] = sra_id_line[1]
                pending_downloads.add(sra_id_line[0])
            except StopIteration:
                download_done = True  # nothing else to download
        build_jobs_int = int(build_jobs[0])
        clean_jobs_int = int(clean_jobs[0])
        if build_jobs_int == 0:  # if we have a build job, we're done - the machine is full
            if to_build_sras and clean_jobs_int == 0 and should_build():
                sra_id, directory = to_build_sras.popitem()
                response['build'] = {'id': sra_id, 'location': directory}
                pending_builds.add(sra_id)
            elif to_clean_sras and (0 <= clean_jobs_int < 4):  # there is room for more clean jobs
                sra_id, dbg_file = to_clean_sras.popitem()
                response['clean'] = {'id': sra_id, 'location': dbg_file}
                pending_cleans.add(sra_id)
        self.send_reply(200, json.dumps(response), {'Content-type': 'application/json'})

    def handle_get_status(self):
        self.send_reply(200, status_str % (
            pending_downloads, pending_builds, pending_cleans, pending_transfers, to_build_sras, to_clean_sras,
            to_transfer_sras, downloaded_sras, built_sras, cleaned_sras, transferred_sras),
                        {'Content-type': 'text/html'})

    def handle_ack(self, operation, post_vars, add_maps, pending_operations):
        sra_id = get_var(post_vars, 'id')
        location = get_var(post_vars, 'location')
        if None in (sra_id, location):
            self.send_reply(400, 'id or location not specified')
            return False
        filename = os.path.join(args.output_dir, f'succeed_{operation}.id')
        with open(filename, 'a') as fp:
            fp.write(f'{sra_id} {location}\n')
        for m in add_maps:
            m[sra_id] = location
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
        self.handle_ack('download', post_vars, [downloaded_sras, to_build_sras], pending_downloads)

    def handle_ack_build(self, post_vars):
        self.handle_ack('build', post_vars, [built_sras, to_clean_sras], pending_builds)

    def handle_ack_clean(self, post_vars):
        if self.handle_ack('clean', post_vars, [cleaned_sras, to_transfer_sras], pending_cleans):
            # a transfer is automatically started after a clean operation by the client
            pending_transfers.add(get_var(post_vars, 'id'))

    def handle_ack_transfer(self, post_vars):
        self.handle_ack('transfer', post_vars, [transferred_sras], pending_transfers)
        del to_transfer_sras[get_var(post_vars, 'id')]

    def handle_nack_download(self, post_vars):
        self.handle_nack('download', post_vars, pending_downloads)

    def handle_nack_build(self, post_vars):
        self.handle_nack('build', post_vars, pending_builds)

    def handle_nack_clean(self, post_vars):
        self.handle_nack('clean', post_vars, pending_cleans)

    def handle_nack_transfer(self, post_vars):
        self.handle_nack('transfer', post_vars, pending_transfers)

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
            self.send_reply(400, 'Bad content-type, only application/x-www-form-urlencoded accepted')
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
        post_vars = self.get_post_vars()
        logging.debug(f'POST {parsed_url} {post_vars}')
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
        else:
            self.send_reply(404, f'Invalid path: {self.path}\n')


def load_file_dict(filename):
    result = {}
    try:
        with open(filename) as fp:
            for line in fp:
                tokens = line.split()
                print(f'{filename} {tokens}')
                assert len(tokens) == 2, 'Invalid line in downloads sras ' + line
                result[tokens[0]] = tokens[1]
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
    global downloaded_sras, built_sras, cleaned_sras, transferred_sras, to_build_sras, to_clean_sras, to_transfer_sras
    downloaded_sras = load_file_dict(filename)
    built_sras = load_file_dict(os.path.join(args.output_dir, 'succeed_build.id'))
    cleaned_sras = load_file_dict(os.path.join(args.output_dir, 'succeed_clean.id'))
    transferred_sras = load_file_dict(os.path.join(args.output_dir, 'succeed_transfer.id'))

    to_build_sras = {k: downloaded_sras[k] for k in set(downloaded_sras) - set(built_sras)}
    to_clean_sras = {k: built_sras[k] for k in set(built_sras) - set(cleaned_sras)}
    to_transfer_sras = {k: cleaned_sras[k] for k in set(cleaned_sras) - set(transferred_sras)}

    # TODO: this incorrect because we need the build paths, not the clean path; need to implement transfer jobs
    to_clean_sras.update(to_transfer_sras)  # because we don't support transfer only jobs
    to_transfer_sras = {}

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
