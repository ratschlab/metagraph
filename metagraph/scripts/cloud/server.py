import argparse
import cgi
import errno
import fileinput
import http.server
import json
import logging
import os
import urllib

# TODO: - maybe handle retries for pending jobs

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


def convert_bucket():
    """ Converts the original gs://... bucket urls in bucket1..bucket9 to sra_id bucket_no pairs for faster parsing """
    logging.info('Converting bucket files to a faster to parse format...')
    fileout = os.path.join(args.data_dir, 'buckets')

    for i in range(1, 9):
        fileout = os.path.join(args.data_dir, f'bucket_proc{i}')
        with open(fileout, 'w') as fpout:
            file = os.path.join(args.data_dir, f'bucket{i}')
            try:
                logging.debug(f'Converting {file}')
                with open(file) as fp:
                    for line in fp:
                        parsed = urllib.parse.urlparse(line)
                        sra_id = parsed.path[1:-2]
                        fpout.write(f'{sra_id}\n')
            except FileNotFoundError:
                logging.fatal(f'Could not find {file} in {args.data_dir}. Please copy it there, otherwise I can\'t '
                              f'figure out which GCloud bucket each sra id is in. Look, I\'m trying to be helpful.')
                exit(1)


class Sra:
    """ Reads a list of files containing SRA ids and returns the ids one by one for further processing """

    def __init__(self, data_files):
        self.data_files = data_files
        self.sraid_to_bucket = {}
        self.add_buckets()

    def next_item(self):
        """Lazy function (generator) to read a sequence of files line by line."""
        for file in self.data_files:
            with open(file) as fp:
                for line in fp:
                    line = line.rstrip().split()
                    if line[0] in downloaded_sras:  # already processed
                        print(f'{line[0]} already downloaded - assuming it\'s ready for processing')
                        continue
                    yield line
        return None

    def load_sraid_to_bucket(self):
        if self.sraid_to_bucket:
            return  # already loaded
        for i in range(1, 9):
            bucket_file = os.path.join(args.data_dir, f'bucket_proc{i}')
            if not os.path.exists(bucket_file):
                convert_bucket()
            logging.debug(f'Loading {bucket_file} ...')
            with open(bucket_file) as fp:
                for line in fp:
                    self.sraid_to_bucket[line.strip()] = i

    def add_buckets(self):
        """ NCBI data on Google cloud is stored in 9 buckets, named gs://sra-pub-run-1..9
            So for each SRA id, we need to know in which bucket to find it. This function looks at the files in
            args.data_dir and if the ids don't contain a bucket id, it adds it.
        """
        if not args.add_gcloud_bucket:
            return
        logging.info('Checking if input files contain bucket information...')
        must_convert = False
        for file in self.data_files:
            with open(file) as fp:
                for line in fp:
                    line = line.rstrip()
                    tokens = line.split()
                    if len(tokens) == 1:
                        must_convert = True
                        break
            if must_convert:
                break
        if not must_convert:
            logging.info('Yup, all good, bucket information is all there')
            return
        logging.info('Bucket information is missing in at least one line. Adding it now...')
        self.load_sraid_to_bucket()
        logging.info('Loaded bucket map')
        total_sras = 0
        skipped_sras = 0
        for file in self.data_files:
            logging.debug(f'Processing {file}')
            for line in fileinput.input(file, inplace=True):
                tokens = line.rstrip().split()
                if len(tokens) == 0:
                    continue
                total_sras += 1
                if not (tokens[0] in self.sraid_to_bucket):
                    logging.warning(f'Sra {tokens[0]} is not found in any bucket, skipping')
                    skipped_sras += 1
                    continue
                print(f'{tokens[0]} {self.sraid_to_bucket[tokens[0]]}')

        logging.info(f'Bucket information added to {len(self.data_files)} files. Skipped {skipped_sras} out of'
                     f' {total_sras} sras.')
        del self.sraid_to_bucket


def get_var(post_vars, var):
    var_list = post_vars.get(var.encode('ascii'))
    if not var_list:
        return None
    return var_list[0].decode('utf-8')


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
                sra_id_line = next(sra_download_gen)
                response['download'] = {'id': sra_id_line[0]}
                if len(sra_id_line) == 2:
                    response['download']['bucket'] = sra_id_line[1]
                pending_downloads.add(sra_id_line[0])
            except StopIteration:
                download_done = True  # nothing else to download
        create_jobs_int = int(create_jobs[0])
        clean_jobs_int = int(clean_jobs[0])
        if create_jobs_int == 0:
            if to_create_sras and clean_jobs_int == 0:
                sra_id, directory = to_create_sras.popitem()
                response['create'] = {'id': sra_id, 'location': directory}
                pending_creates.add(sra_id)
            elif to_clean_sras and (0 <= clean_jobs_int < 4):  # there is room for more clean jobs
                sra_id, dbg_file = to_clean_sras.popitem()
                response['clean'] = {'id': sra_id, 'location': dbg_file}
                pending_cleans.add(sra_id)
        self.send_reply(200, json.dumps(response), {'Content-type': 'application/json'})

    def handle_get_status(self):
        self.send_reply(200, status_str % (
            pending_downloads, pending_creates, pending_cleans, pending_transfers, to_create_sras, to_clean_sras,
            to_transfer_sras, downloaded_sras, created_sras, cleaned_sras, transferred_sras),
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
        self.handle_ack('download', post_vars, [downloaded_sras, to_create_sras], pending_downloads)

    def handle_ack_create(self, post_vars):
        self.handle_ack('create', post_vars, [created_sras, to_clean_sras], pending_creates)

    def handle_ack_clean(self, post_vars):
        if self.handle_ack('clean', post_vars, [cleaned_sras, to_transfer_sras], pending_cleans):
            # a transfer is automatically started after a clean operation by the client
            pending_transfers.add(get_var(post_vars, 'id'))

    def handle_ack_transfer(self, post_vars):
        self.handle_ack('transfer', post_vars, [transferred_sras], pending_transfers)
        del to_transfer_sras[get_var(post_vars, 'id')]

    def handle_nack_download(self, post_vars):
        self.handle_nack('download', post_vars, pending_downloads)

    def handle_nack_create(self, post_vars):
        self.handle_nack('create', post_vars, pending_creates)

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
            self.send_reply(400, "Bad content-type, only application/x-www-form-urlencoded accepted")
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
        if not post_vars:
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
    global sra_download_gen
    sra_download_gen = Sra(data_files).next_item()

    filename = os.path.join(args.output_dir, 'succeed_download.id')
    if not os.path.exists(os.path.dirname(filename)):
        try:
            os.makedirs(os.path.dirname(filename))
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    global downloaded_sras, created_sras, cleaned_sras, transferred_sras, to_create_sras, to_clean_sras, to_transfer_sras
    downloaded_sras = load_file_dict(filename)
    created_sras = load_file_dict(os.path.join(args.output_dir, 'succeed_create.id'))
    cleaned_sras = load_file_dict(os.path.join(args.output_dir, 'succeed_clean.id'))
    transferred_sras = load_file_dict(os.path.join(args.output_dir, 'succeed_transfer.id'))

    to_create_sras = {k: downloaded_sras[k] for k in set(downloaded_sras) - set(created_sras)}
    to_clean_sras = {k: created_sras[k] for k in set(created_sras) - set(cleaned_sras)}
    to_transfer_sras = {k: cleaned_sras[k] for k in set(cleaned_sras) - set(transferred_sras)}

    # TODO: this incorrect because we need the create paths, not the clean path; need to implement transfer jobs
    to_clean_sras.update(to_transfer_sras)  # because we don't support transfer only jobs
    to_transfer_sras = {}


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

    global args
    args = parser.parse_args()


if __name__ == '__main__':
    parse_args()

    init_logging()

    init_state()

    httpd = http.server.HTTPServer(('', args.port), SimpleHTTPRequestHandler)
    logging.info(f'Starting server on port {args.port}')
    httpd.serve_forever()
