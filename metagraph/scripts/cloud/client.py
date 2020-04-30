#!/usr/bin/env python3

import argparse
import http
import json
import logging
import os
import socket
import subprocess
import time
import urllib.request

args = None

download_processes = {}
build_processes = {}
clean_processes = {}
transfer_processes = {}

# the location (host:/path) where clean and build get their files from; used to clean up files that are processed
build_source = {}
clean_source = {}

MAX_DOWNLOAD_PROCESSES = 1
MAX_BUILD_PROCESSES = 1
MAX_CLEAN_PROCESSES = 4


def get_work():
    url = f'http://{args.server_host}:{args.server_port}/jobs?download={len(download_processes)}&' \
          f'build={len(build_processes)}&clean={len(clean_processes)}&client_id={args.client_id}'
    for i in range(10):
        try:
            response = urllib.request.urlopen(url)
            if response.getcode() == 200:
                return json.loads(response.read())
            elif response.getcode() == 204:
                return None
            else:
                logging.warning(f'Server returned response code {response.getcode()} for {url}')
                time.sleep(5)  # avoid overwhelming the server
        except urllib.error.URLError as e:
            logging.error(f'Failed to open URL {url} Reason: {e.reason}')
            time.sleep(5)  # wait a bit and try again
        except http.client.RemoteDisconnected as e:
            logging.error(f'Failed to open URL {url} Reason: remote disconnected.')
            time.sleep(5)  # wait a bit and try again
    return {}


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


def clean_dir():
    return os.path.join(args.output_dir, 'cleaned_graphs/')


def clean_file(sra_id):
    return os.path.join(clean_dir(), f'{sra_id}.*')


def start_download(download_resp):
    sra_id = download_resp['id']
    if args.source == 'ena':
        download_processes[sra_id] = subprocess.Popen(['./download_ena.sh', sra_id, download_dir_base()])
    else:
        download_processes[sra_id] = subprocess.Popen(
            ['./download.sh', download_resp['bucket'], sra_id, download_dir_base()])


def internal_ip():
    try:
        return socket.gethostbyname(socket.gethostname())
    except socket.gaierror:
        return '127.0.0.1'  # this usually happens on dev laptops; cloud machines work fine


def start_build(sra_id, location):
    input_dir_base = download_dir_base()
    if not location.startswith(internal_ip()):
        build_source[sra_id] = location
        return_code = subprocess.call(['scp', '-r', location, input_dir_base])
        if return_code != 0:
            logging.warning(f'Copying from {location} to {input_dir_base} failed')
            return False
        else:
            logging.info(f'Copying from {location} to {input_dir_base} completed successfully')
    else:
        logging.info(f'Luckily {location} already is on this machine.')
    input_dir = download_dir(sra_id)
    output_dir = build_dir(sra_id)
    build_processes[sra_id] = subprocess.Popen(['./build.sh', sra_id, input_dir, output_dir])
    return True


def make_dir_if_needed(path):
    try:
        os.makedirs(path, exist_ok=True)
    except FileExistsError:
        pass


def start_clean(sra_id, location):
    input_dir = build_dir_base()
    if not location.startswith(internal_ip()):
        clean_source[sra_id] = location
        return_code = subprocess.call(['scp', '-r', location, input_dir])
        if return_code != 0:
            logging.warning(f'Copying from {location} to {input_dir} failed')
            return False
        else:
            logging.info(f'Copying from {location} to {input_dir} completed successfully')
    else:
        logging.info('Luckily the files to use for building the graph already are on this machine.')
    input_file = build_file(sra_id)
    output_file = os.path.join(clean_dir(), sra_id)
    clean_processes[sra_id] = subprocess.Popen(['./clean.sh', sra_id, input_file, output_file])
    return True


def start_transfer(sra_id, cleaned_graph_location):
    transfer_processes[sra_id] = subprocess.Popen(
        f'scp -r {cleaned_graph_location} {args.destination}; rm -rf {cleaned_graph_location}', shell=True)


def is_full():
    return len(download_processes) == MAX_DOWNLOAD_PROCESSES and (
            len(clean_processes) == MAX_CLEAN_PROCESSES or len(build_processes) == MAX_BUILD_PROCESSES)


def ack(operation, sra_id, location):
    url = f'http://{args.server_host}:{args.server_port}/jobs/ack/{operation}'
    data = f'id={sra_id}&location={location}&client_id={args.client_id}'
    while True:
        try:
            request = urllib.request.Request(url, data=data.encode('UTF-8'),
                                             headers={'Content-type': 'application/x-www-form-urlencoded'},
                                             method='POST')
            response = urllib.request.urlopen(request)
            if response.getcode() == 200:
                return True
            elif response.getcode() == 204:
                return False
            else:
                logging.warning(f'Server returned response code {response.getcode()} for {url}')
                time.sleep(5)  # avoid overwhelming the server
        except urllib.error.URLError as e:
            logging.error(f'Failed to open URL {url} Reason: {e.reason}')
            time.sleep(5)  # wait a bit and try again


def nack(operation, sra_id):
    url = f'http://{args.server_host}:{args.server_port}/jobs/nack/{operation}'
    data = f'id={sra_id}&client_id={args.client_id}'
    while True:
        try:
            request = urllib.request.Request(url, data=data.encode('UTF-8'),
                                             headers={'Content-type': 'application/x-www-form-urlencoded'},
                                             method='POST')
            response = urllib.request.urlopen(request)
            if response.getcode() == 200:
                return True
            elif response.getcode() == 204:
                return False
            else:
                logging.warning(f'Server returned response code {response.getcode()} for {url}')
                time.sleep(5)  # avoid overwhelming the server
        except urllib.error.URLError as e:
            logging.error(f'Failed to open URL {url} Reason: {e.reason}')
            time.sleep(5)  # wait a bit and try again


def check_status():
    completed_downloads = set()
    for sra_id, download_process in download_processes.items():
        return_code = download_process.poll()
        if return_code is not None:
            completed_downloads.add(sra_id)
            if return_code == 0:
                logging.info(f'Download for SRA id {sra_id} completed successfully.')
                location = os.path.join(f'{internal_ip()}:{download_dir_base()}', sra_id)

                if not ack('download', sra_id, location):  # all done, yay!
                    return False;
            else:
                logging.warning(f'Download for SRA id {sra_id} failed.')
                if not nack('download', sra_id):
                    return False  # this shouldn't happen, as we still have work to do
    for d in completed_downloads:
        del download_processes[d]

    completed_builds = set()
    for sra_id, build_process in build_processes.items():
        return_code = build_process.poll()

        if return_code is not None:
            completed_builds.add(sra_id)
            # clean up the download path; if adding retries, do this only on success
            if sra_id in build_source:
                download_path = build_source[sra_id]
                logging.info(f'Cleaning up {download_path}')
                parsed_download_path = urllib.parse.urlparse(download_path)
                subprocess.run(['ssh', parsed_download_path.scheme, f'rm -rf {parsed_download_path.path}'])
                del build_source[sra_id]

            if return_code == 0:
                logging.info(f'Building graph for SRA id {sra_id} completed successfully.')
                location = f'{internal_ip()}:{build_dir(sra_id)}'

                if not ack('build', sra_id, location):  # all done, yay!
                    return False;
            else:
                logging.warning(f'Building graph for SRA id {sra_id} failed.')
                if not nack('build', sra_id):
                    return False  # this shouldn't happen, as we still have work to do
    for d in completed_builds:
        del build_processes[d]

    completed_cleans = set()
    for sra_id, clean_process in clean_processes.items():
        return_code = clean_process.poll()
        if return_code is not None:
            completed_cleans.add(sra_id)

            if sra_id in clean_source:
                # clean up the original graph; if adding retries, do this only on success
                clean_path = clean_source[sra_id]
                logging.info(f'Cleaning up {clean_path}')
                parsed_download_path = urllib.parse.urlparse(clean_path)
                subprocess.run(['ssh', parsed_download_path.scheme, f'rm -rf {parsed_download_path.path}'])
                del clean_source[sra_id]

            if return_code == 0:
                logging.info(f'Cleaning graph for SRA id {sra_id} completed successfully.')
                location = f'{internal_ip()}:{clean_file(sra_id)}'

                if not ack('clean', sra_id, location):  # all done, yay!
                    return False


                start_transfer(sra_id, clean_file(sra_id))
            else:
                logging.warning(f'Cleaning graph for SRA id {sra_id} failed.')
                if not nack('clean', sra_id):
                    return False  # this shouldn't happen, as we still have work to do
    for d in completed_cleans:
        del clean_processes[d]

    completed_transfers = set()
    for sra_id, transfer_process in transfer_processes.items():
        return_code = transfer_process.poll()
        if return_code is not None:
            completed_transfers.add(sra_id)
            if return_code == 0:
                logging.info(f'Transferring graph for SRA id {sra_id} completed successfully.')
                if not ack('transfer', sra_id, args.destination):  # all done, yay!
                    return False
            else:
                logging.warning(f'Transferring cleaned graph for SRA id {sra_id} failed.')
                if not nack('transfer', sra_id):
                    return False  # this shouldn't happen, as we still have work to do
    for d in completed_transfers:
        del transfer_processes[d]

    return True


def do_work():
    while True:
        if not check_status():
            break
        if is_full():
            time.sleep(5)
            continue
        work_response = get_work()
        logging.debug(f'Response is {work_response}')
        if work_response is None:
            break
        if 'download' in work_response:
            start_download(work_response['download'])
        if 'build' in work_response:
            if not start_build(work_response['build']['id'], work_response['build']['location']):
                nack('build', work_response['build']['id'])
        if 'clean' in work_response:
            if not start_clean(work_response['clean']['id'], work_response['clean']['location']):
                nack('clean', work_response['clean']['id'])
        time.sleep(5)


def check_env():
    """ Make sure all the necessary software is in place to successfully run the client and create working
    directories """

    make_dir_if_needed(download_dir_base())
    make_dir_if_needed(build_dir_base())
    make_dir_if_needed(clean_dir())

    logging.basicConfig(format='%(asctime)s %(message)s', level=logging.DEBUG)
    file_handler = logging.FileHandler("{0}/{1}.log".format(args.output_dir, 'client'))
    file_handler.setLevel(logging.DEBUG)
    logging.getLogger().addHandler(file_handler)

    if subprocess.call(['./prereq.sh']) != 0:
        logging.error("Some prerequisites are missing on this machine. Bailing out.")
        exit(1)

    make_dir_if_needed(download_dir_base())
    make_dir_if_needed(build_dir_base())
    make_dir_if_needed(clean_dir())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--source', help='Where to download the data from: ena or ncbi', choices=('ena', 'ncbi'),
                        default='ncbi')
    parser.add_argument('--server_host', help='HTTP server name or ip')
    parser.add_argument('--server_port', default=8000, help='HTTP Port on which the server runs')
    parser.add_argument(
        '--output_dir',
        default=os.path.expanduser('~/.metagraph/'),
        help='Location of the directory containing the input data')
    parser.add_argument('--destination', default='leomed:/cluster/work/grlab/projects/metagenome/scratch/cloud',
                        help='Host/directory where the cleaned BOSS graphs are copied to')
    parser.add_argument('--client_id', default='-1',
                        help='Unique id for each client, used to identify clients for logging purposes')
    args = parser.parse_args()

    if not os.path.isabs(args.output_dir):
        logging.error(f'output_dir must be an absolute path, not {args.output_dir}')
        exit(1)
    if not args.server_host:
        logging.error('missing --server_host. Can\'t connect to server without it!')
        exit(1)
    check_env()

    do_work()
