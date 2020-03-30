#!/usr/bin/env python3
# TODO: get a valid client id, perhaps from managed cluster

import argparse
import collections
import http
import http.server
import json
import logging
import os
import psutil
import signal
import socket
import subprocess
import threading
import time
import urllib.parse
import urllib.request

args = None

sra_info = {}  # global information about a processed SRA (for now, time only)

download_processes = {}
build_processes = {}
clean_processes = {}
transfer_processes = {}

waiting_builds = collections.OrderedDict()
waiting_cleans = collections.OrderedDict()

MAX_DOWNLOAD_PROCESSES = 1
MAX_build_PROCESSES = 1
MAX_CLEAN_PROCESSES = 4

downloads_done = False
must_quit = False

status_str = f"""
<html>
<head>
<title>Status of metagraph server</title>
</head>
<body>
<h3> Currently running jobs </h3>
<p>Downloading: %s</p>
<p>Creating: %s</p>
<p>Cleaning: %s</p>
<p>Transferring: %s</p>

<h3> Waiting </h3>
<p>Waiting builds %s</p>
<p>Waiting cleans %s</p>

<h3> Download status </h3>
Download done: %s

</body>
</html>
"""


def get_work():
    global downloads_done
    if downloads_done or len(download_processes) > 0 or len(waiting_builds) > 1:
        return None
    url = f'http://{args.server}/jobs'
    for i in range(10):
        if must_quit:
            break
        try:
            response = urllib.request.urlopen(url, timeout=15)
            if response.getcode() == 200:
                response = json.loads(response.read())
                logging.debug(f'Response is: {response}')
                if 'download' not in response:
                    downloads_done = True
                    return None
                return response
            else:
                logging.warning(f'Server returned response code {response.getcode()} for {url}')
                time.sleep(5)  # avoid overwhelming the server
        except urllib.error.URLError as e:
            logging.error(f'Failed to open URL {url} Reason: {e.reason}')
            time.sleep(5)  # wait a bit and try again
        except http.client.RemoteDisconnected as e:
            logging.error(f'Failed to open URL {url} Reason: remote disconnected.')
            time.sleep(5)  # wait a bit and try again
        except:
            logging.error(f'Failed to open URL {url} Reason: unknown.')
            time.sleep(5)  # wait a bit and try again

    return None


def download_dir_base():
    return os.path.join(args.output_dir, 'downloads/')


def download_dir(sra_id):
    return os.path.join(download_dir_base(), sra_id)


def build_dir_base():
    return os.path.join(args.output_dir, 'graphs')


def build_dir(sra_id):
    return os.path.join(build_dir_base(), sra_id)


def build_file(sra_id):
    return os.path.join(build_dir(sra_id), f'{sra_id}.dbg')


def clean_dir_base():
    return os.path.join(args.output_dir, 'cleaned_graphs')


def clean_dir(sra_id):
    return os.path.join(clean_dir_base(), sra_id)


def start_download(download_resp):
    if 'id' not in download_resp:
        logging.info('No more downloads available. We\'re almost done!')
        return
    sra_id = download_resp['id']
    if args.source == 'ena':
        download_processes[sra_id] = (subprocess.Popen(['./download_ena.sh', sra_id, download_dir_base()]), time.time())
    else:
        if not download_resp['bucket']:
            logging.fatal('Specified NCBI as download source, but server response has no "bucket" field')
        make_dir_if_needed(download_dir(sra_id))
        log_file_name = os.path.join(download_dir(sra_id), 'download.log')
        bucket = download_resp['bucket']
        download_processes[sra_id] = (subprocess.Popen(
            f'./download_ncbi.sh {bucket} {sra_id} {download_dir_base()}  2>&1 | '
            f'grep -v "Stage" >{log_file_name}', shell=True),
                                      time.time())
    sra_info[sra_id] = (time.time(),)


def internal_ip():
    try:
        return socket.gethostbyname(socket.gethostname())
    except socket.gaierror:
        return '127.0.0.1'  # this usually happens on dev laptops; cloud machines work fine


def start_build(sra_id, wait_time, buffer_size_gb):
    input_dir = download_dir(sra_id)
    output_dir = build_dir(sra_id)
    logging.info(f'Starting build from {input_dir} to {output_dir}, ram={round(buffer_size_gb, 2)}')
    make_dir_if_needed(build_dir(sra_id))
    log_file_name = os.path.join(build_dir(sra_id), 'build.log')
    build_processes[sra_id] = (subprocess.Popen(
        f'./build.sh {sra_id} {input_dir} {output_dir} {buffer_size_gb} 2>&1 > {log_file_name}',
        shell=True), time.time(), wait_time)
    return True


def make_dir_if_needed(path):
    try:
        os.makedirs(path, exist_ok=True)
    except FileExistsError:
        pass


def start_clean(sra_id, wait_time):
    input_file = build_file(sra_id)
    make_dir_if_needed(clean_dir(sra_id))
    output_file = os.path.join(clean_dir(sra_id), sra_id)
    logging.info(f'Starting clean with {input_file} {output_file}')
    clean_file_name = os.path.join(clean_dir(sra_id), 'clean.log')
    clean_processes[sra_id] = (
        subprocess.Popen(f'./clean.sh {sra_id} {input_file} {output_file} 2>&1 | grep -v "%," > {clean_file_name}',
                         shell=True),
        time.time(), wait_time)
    return True


def start_transfer(sra_id, cleaned_graph_location):
    transfer_processes[sra_id] = (subprocess.Popen(
        f'gsutil -q -u metagraph cp -r {cleaned_graph_location} {args.destination}', shell=True), time.time())


def ack(operation, params):
    url = f'http://{args.server}/jobs/ack/{operation}'
    data = f'client_id={internal_ip()}&' + urllib.parse.urlencode(params)
    while True:
        try:
            request = urllib.request.Request(url, data=data.encode('UTF-8'),
                                             headers={'Content-type': 'application/x-www-form-urlencoded'},
                                             method='POST')
            response = urllib.request.urlopen(request)
            if response.getcode() == 200:
                return
            else:
                logging.warning(f'Server returned response code {response.getcode()} for {url}')
                time.sleep(5)  # avoid overwhelming the server
        except urllib.error.URLError as e:
            logging.error(f'Failed to open URL {url} Reason: {e.reason}')
            time.sleep(5)  # wait a bit and try again
        except http.client.RemoteDisconnected as e:
            logging.error(f'Failed to open URL {url} Reason: remote disconnected.')
            time.sleep(5)  # wait a bit and try again
        except:
            logging.error(f'Failed to open URL {url} Reason: unknown.')
            time.sleep(5)  # wait a bit and try again


def nack(operation, params):
    url = f'http://{args.server}/jobs/nack/{operation}'
    data = f'client_id={internal_ip()}&' + urllib.parse.urlencode(params)
    while True:
        try:
            request = urllib.request.Request(url, data=data.encode('UTF-8'),
                                             headers={'Content-type': 'application/x-www-form-urlencoded'},
                                             method='POST')
            response = urllib.request.urlopen(request)
            if response.getcode() == 200:
                return
            else:
                logging.warning(f'Server returned response code {response.getcode()} for {url}')
                time.sleep(5)  # avoid overwhelming the server
        except urllib.error.URLError as e:
            logging.error(f'Failed to open URL {url} Reason: {e.reason}')
            time.sleep(5)  # wait a bit and try again
        except http.client.RemoteDisconnected as e:
            logging.error(f'Failed to open URL {url} Reason: remote disconnected.')
            time.sleep(5)  # wait a bit and try again
        except:
            logging.error(f'Failed to open URL {url} Reason: unknown.')
            time.sleep(5)  # wait a bit and try again


def dir_size(dir_path):
    total_size = 0
    for dir_path, dir_names, file_names in os.walk(dir_path):
        for f in file_names:
            fp = os.path.join(dir_path, f)
            # skip if it is symbolic link
            if not os.path.islink(fp):
                total_size += os.path.getsize(fp)

    return round(total_size / 1000000, 2)


def check_status():
    global must_quit
    if must_quit:
        return False
    completed_downloads = set()
    for sra_id, (download_process, start_time) in download_processes.items():
        return_code = download_process.poll()
        is_timed_out = (time.time() - start_time) > 15 * 60
        if return_code is not None or is_timed_out:
            completed_downloads.add(sra_id)
            log_file_name = os.path.join(download_dir(sra_id), 'download.log')
            logging.info(f'Download finished with output\n {open(log_file_name).read()}')
            download_path = download_dir(sra_id)
            sra_dir = os.path.join(download_path, 'sra')
            download_size_mb = dir_size(sra_dir)
            subprocess.run(['rm', '-rf', sra_dir])
            kmc_dir = os.path.join(download_path, 'kmc')
            kmc_size_mb = dir_size(kmc_dir)
            success = True
            if return_code == 0 or return_code == 3:
                if return_code == 0:
                    logging.info(f'Download for SRA id {sra_id} completed successfully.')
                else:
                    logging.info(f'Download for SRA id {sra_id} completed with some errors.')
                stats_file = os.path.join(download_path, 'stats')
                try:
                    with open(stats_file) as stats:
                        json_resp = json.loads(stats.read())
                    if 'Stats' in json_resp and '#Unique_counted_k-mers' in json_resp['Stats']:
                        kmer_count_unique = int(json_resp['Stats']['#Unique_counted_k-mers'])
                        kmer_count_total = int(json_resp['Stats']['#Total no. of k-mers'])
                    else:
                        logging.warning('KMC returned no stats, assuming failure')
                        success = False
                except FileNotFoundError:
                    logging.warning(f'Could not find KMC stats file {stats_file}, baling out.')
                    success = False
            else:
                success = False
            if success:
                params = {'id': sra_id, 'time': int(time.time() - start_time), 'size_mb': download_size_mb,
                          'kmc_size_mb': kmc_size_mb, 'kmer_count_unique': kmer_count_unique,
                          'kmer_count_total': kmer_count_total}
                sra_info[sra_id] = (*sra_info[sra_id], download_size_mb, kmer_count_unique)
                ack('download', params)
                waiting_builds[sra_id] = (time.time())
            else:
                logging.warning(f'Download for SRA id {sra_id} failed. Removing {download_path}')
                subprocess.run(['rm', '-rf', download_path])
                params = {'id': sra_id, 'time': int(time.time() - start_time), 'size_mb': download_size_mb,
                          'kmc_size_mb': kmc_size_mb}
                nack('download', params)
    for d in completed_downloads:
        del download_processes[d]

    completed_builds = set()
    for sra_id, (build_process, start_time, wait_time) in build_processes.items():
        return_code = build_process.poll()

        if return_code is not None:
            completed_builds.add(sra_id)
            log_file_name = os.path.join(build_dir(sra_id), 'build.log')
            logging.info(f'Build finished with output\n {open(log_file_name).read()}')

            # clean up the download path; if adding retries, do this only on success
            download_path = download_dir(sra_id)
            logging.info(f'Cleaning up {download_path}')
            subprocess.run(['rm', '-rf', download_path])

            build_path = build_dir(sra_id)
            build_size = dir_size(build_path)
            if return_code == 0:
                logging.info(f'Building graph for SRA id {sra_id} completed successfully.')
                params = {'id': sra_id, 'time': int(time.time() - start_time),
                          'wait_time': int(wait_time), 'size_mb': build_size}
                ack('build', params)
                waiting_cleans[sra_id] = (time.time())
            else:
                logging.warning(f'Building graph for SRA id {sra_id} failed. Removing {build_path}.')
                subprocess.run(['rm', '-rf', build_path])
                params = {'id': sra_id, 'time': int(time.time() - start_time), 'wait_time': int(wait_time),
                          'size_mb': build_size, 'return_code': return_code}
                nack('build', params)
    for d in completed_builds:
        del build_processes[d]

    completed_cleans = set()
    for sra_id, (clean_process, start_time, wait_time) in clean_processes.items():
        return_code = clean_process.poll()
        if return_code is not None:
            completed_cleans.add(sra_id)
            log_file_name = os.path.join(clean_dir(sra_id), 'clean.log')
            logging.info(f'Clean finished with output\n {open(log_file_name).read()}')

            # clean up the original graph; if adding retries, do this only on success
            build_path = build_dir(sra_id)
            logging.info(f'Cleaning up {build_path}')
            subprocess.run(['rm', '-rf', build_path])

            cleaned_dir = clean_dir(sra_id)
            cleaned_size = dir_size(cleaned_dir)
            if return_code == 0:
                logging.info(f'Cleaning graph for SRA id {sra_id} completed successfully.')

                params = {'id': sra_id, 'time': int(time.time() - start_time),
                          'size_mb': cleaned_size}
                ack('clean', params)
                start_transfer(sra_id, cleaned_dir)
            else:
                logging.warning(f'Cleaning graph for SRA id {sra_id} failed. Removing {cleaned_dir}')
                subprocess.run(['rm', '-rf', cleaned_dir])
                params = {'id': sra_id, 'time': int(time.time() - start_time),
                          'size_mb': cleaned_size, 'return_code': return_code}
                nack('clean', params)
    for d in completed_cleans:
        del clean_processes[d]

    completed_transfers = set()
    for sra_id, (transfer_process, start_time) in transfer_processes.items():
        return_code = transfer_process.poll()
        if return_code is not None:
            completed_transfers.add(sra_id)
            # clean up the cleaned graph; if adding retries, do this only on success
            clean_path = clean_dir(sra_id)
            cleaned_size = dir_size(clean_path)
            logging.info(f'Cleaning up {clean_path}')
            subprocess.run(['rm', '-rf', clean_path])

            if return_code == 0:
                logging.info(f'Transferring graph for SRA id {sra_id} completed successfully.')
                params = {'id': sra_id, 'time': int(time.time() - start_time),
                          'total_time': int(time.time() - sra_info[sra_id][0]), 'size_init_mb': sra_info[sra_id][1],
                          'size_final_mb': cleaned_size}
                ack('transfer', params)
            else:
                logging.warning(f'Transferring cleaned graph for SRA id {sra_id} failed.')
                params = {'id': sra_id, 'time': int(time.time() - start_time),
                          'size_mb': cleaned_size}
                nack('transfer', params)

    available_ram_gb = psutil.virtual_memory().available / 1e9 - 1
    if len(build_processes) == 0 and len(clean_processes) < 2 and waiting_builds and len(waiting_cleans) < 3:
        sra_id, (start_time) = waiting_builds.popitem()
        num_kmers = sra_info[sra_id][2]
        required_ram_gb = round((num_kmers * 2) / 1e9 + 1, 2)
        total_ram_gb = psutil.virtual_memory().total / 1e9
        if required_ram_gb > total_ram_gb - 3:
            build_path = build_dir(sra_id)
            logging.warning(f'Building graph for SRA id {sra_id} needs too much RAM '
                            f'({required_ram_gb}GB). Removing {build_path}.')
            subprocess.run(['rm', '-rf', build_path])
            params = {'id': sra_id, 'time': int(time.time() - start_time), 'wait_time': int(wait_time),
                      'size_mb': build_size, 'required_ram_gb': required_ram_gb}
            nack('build', params)
        elif required_ram_gb < available_ram_gb and available_ram_gb > 2:
            buffer_size_gb = max(2, min(round(required_ram_gb * 0.8 - 1), 20))
            start_build(sra_id, time.time() - start_time, buffer_size_gb)
        else:
            logging.info(f'Not enough RAM for building. Have {round(available_ram_gb, 2)}GB need {required_ram_gb}GB')

    # we do max one clean while we have a build, and up to 4 cleans if no build is running
    if (len(clean_processes) == 0 or (len(clean_processes) < 3 and len(build_processes) == 0)) and waiting_cleans:
        available_ram_gb = psutil.virtual_memory().available / 1e9 - 1
        for sra_id, (start_time) in waiting_cleans.items():
            # remove the old clean waiting and append the new one after
            build_path = build_dir(sra_id)
            build_size_gb = dir_size(build_path) / 1e9
            required_ram_gb = max(build_size_gb * 1.1, build_size_gb + 1)

            if available_ram_gb - required_ram_gb > 0:
                start_clean(sra_id, time.time() - start_time)
                del waiting_cleans[sra_id]
                break
            logging.info(f'Not enough RAM for cleaning {sra_id}. '
                         f'Have {round(available_ram_gb, 2)}GB need {round(build_size_gb + 0.5, 2)}GB')

    for d in completed_transfers:
        del transfer_processes[d]
    return download_processes or build_processes or clean_processes or transfer_processes or not downloads_done


def do_work():
    i = 0
    while True:
        if not check_status():
            break

        work_response = get_work()
        if work_response is None:  # already downloading/waiting to build or no more downloads
            time.sleep(10)
            continue
        if 'download' in work_response:
            start_download(work_response['download'])
        else:
            logging.error(f'Server response invalid. Expected a \'download\' tag: {work_response}')
        if i % 10 == 0:
            subprocess.Popen(
                [f'gsutil rsync -x \'(?!^client.log$)\' {args.output_dir} {args.log_destination}/{internal_ip()}/'],
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)

        i = i + 1
        time.sleep(5)


def check_env():
    """ Make sure all the necessary software is in place to successfully run the client and create working
    directories """

    make_dir_if_needed(download_dir_base())
    make_dir_if_needed(build_dir_base())
    make_dir_if_needed(clean_dir_base())

    logging.basicConfig(format='%(asctime)s %(message)s', level=logging.DEBUG, datefmt='%Y-%m-%d %H:%M:%S')
    file_handler = logging.FileHandler(f'{args.output_dir}/client.log')
    file_handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s %(message)s')
    file_handler.setFormatter(formatter)
    logging.getLogger().addHandler(file_handler)

    if subprocess.call(['./prereq.sh']) != 0:
        logging.error("Some prerequisites are missing on this machine. Bailing out.")
        exit(1)


def handle_quit():
    ids = ','.join(
        list(download_processes) + list(build_processes) + list(clean_processes) + list(transfer_processes) + list(
            waiting_builds) + list(waiting_cleans))
    url = f'http://{args.server}/jobs/preempt'
    data = f'client_id={internal_ip()}&ids={ids}'
    while True:
        try:
            request = urllib.request.Request(url, data=data.encode('UTF-8'),
                                             headers={'Content-type': 'application/x-www-form-urlencoded'},
                                             method='POST')
            response = urllib.request.urlopen(request, timeout=10)
            if response.getcode() == 200:
                break
            else:
                logging.warning(f'Server returned response code {response.getcode()} for {url}')
                time.sleep(5)  # avoid overwhelming the server
        except urllib.error.URLError as e:
            logging.error(f'Failed to open URL {url} Reason: {e.reason}')
            time.sleep(5)  # wait a bit and try again
        except http.client.RemoteDisconnected as e:
            logging.error(f'Failed to open URL {url} Reason: remote disconnected.')
            time.sleep(5)  # wait a bit and try again
        except:
            logging.error(f'Failed to open URL {url} Reason: unknown.')
            time.sleep(5)  # wait a bit and try again
    global must_quit
    must_quit = True
    for k, v in download_processes.items():
        v[0].kill()
    for k, v in build_processes.items():
        v[0].kill()
    for k, v in clean_processes.items():
        v[0].kill()
    for k, v in transfer_processes.items():
        v[0].kill()


class SimpleHTTPRequestHandler(http.server.BaseHTTPRequestHandler):
    """ Processes status and 'please die' requests """

    def handle_get_status(self):
        self.send_reply(200, status_str % (download_processes,
                                           build_processes,
                                           clean_processes,
                                           transfer_processes,
                                           waiting_builds,
                                           waiting_cleans,
                                           downloads_done),
                        {'Content-type': 'text/html'})

    def do_GET(self):
        parsed_url = urllib.parse.urlparse(self.path)
        if parsed_url.path == '/quit':
            print('SRA Client was asked to quit. Notifying server.')
            handle_quit()
            print('Good-bye.')
            exit(0)
        elif parsed_url.path == '/status':
            self.handle_get_status()
        elif parsed_url.path == '/healthz':
            if not must_quit:
                self.send_reply(200, 'OK')
            else:
                self.send_repy(200, 'Quit')
        else:
            self.send_reply(404, f'Invalid path: {self.path}\n')

    def send_reply(self, code, message, headers={}):
        self.send_response(code)
        for k, v in headers.items():
            self.send_header(k, v)
        self.end_headers()
        self.flush_headers()
        self.wfile.write(message.encode('utf-8'))
        self.wfile.flush()


def signal_handler(signum, frame):
    print(f'Received signal {signum}. Gracefully dying.')
    handle_quit()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--source', help='Where to download the data from: ena or ncbi', choices=('ena', 'ncbi'),
                        default='ncbi')
    parser.add_argument('--server', help='HTTP server host:port')
    parser.add_argument(
        '--output_dir',
        default=os.path.expanduser('~/.metagraph/'),
        help='Location of the directory containing the input data')
    parser.add_argument('--destination', default='gs://metagraph14/clean/',
                        help='Host/directory where the cleaned BOSS graphs are copied to')
    parser.add_argument('--log_destination', default='gs://metagraph14/logs',
                        help='GS folder where client logs are collected')
    parser.add_argument('--port', default=8001, help='HTTP Port on which the status/kill server runs')
    args = parser.parse_args()
    if not os.path.isabs(args.output_dir):
        logging.error(f'output_dir must be an absolute path, not {args.output_dir}')
        exit(1)

    check_env()

    if not args.server:
        logging.info('Trying to find server address...')
        if subprocess.call(['gsutil', 'cp', 'gs://metagraph14/server', '/tmp/server'], stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE) != 0:
            logging.error('Cannot find server ip/port on Google Cloud Storage. Sorry, I tried.')
            exit(1)
        with open('/tmp/server') as fp:
            args.server = fp.read()

        logging.info(f'Found server at {args.server}')

    # gracefully handle termination
    signal.signal(signal.SIGTERM, signal_handler)

    httpd = http.server.HTTPServer(('', args.port), SimpleHTTPRequestHandler)
    logging.info(f'Starting client status server on port {args.port}')
    thread = threading.Thread(name='server_thread', target=httpd.serve_forever)
    thread.daemon = True
    thread.start()

    do_work()
