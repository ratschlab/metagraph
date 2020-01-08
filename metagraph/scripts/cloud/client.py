import argparse
import json
import logging
import os
import pathlib
import socket
import subprocess
import time
import urllib.request

logger = logging.getLogger('metagraph-client')
args = None

download_processes = {}
create_processes = {}
clean_processes = {}
transfer_processes = {}

MAX_DOWNLOAD_PROCESSES = 1
MAX_CREATE_PROCESSES = 1
MAX_CLEAN_PROCESSES = 4


def get_work():
    url = f'http://{args.server_host}:{args.server_port}/jobs?download={len(download_processes)}&' \
          f'create={len(create_processes)}&clean={len(clean_processes)}'
    for i in range(10):
        try:
            response = urllib.request.urlopen(url)
            if response.getcode() == 200:
                return json.loads(response.read())
            elif response.getcode() == 204:
                return None
            else:
                logger.warning(f'Server returned response code {response.getcode()} for {url}')
                time.sleep(5)  # avoid overwhelming the server
        except urllib.error.URLError as e:
            print(f'Failed to open URL {url} Reason: {e.reason}')
            time.sleep(5)  # wait a bit and try again
    return {}


def download_dir():
    return f'{args.output_dir}/downloads/'


def create_dir():
    return f'{args.output_dir}/graphs/'


def create_file(sra_id):
    return os.path.join(create_dir(), f'{sra_id}.dbg')


def clean_dir():
    return f'{args.output_dir}/cleaned_graphs/'


def clean_file(sra_id):
    return os.path.join(clean_dir(), f'{sra_id}.dbg')  # TODO: figure out the proper file name


def start_download(sra_id):
    download_processes[sra_id] = subprocess.Popen(['./download.sh', sra_id, download_dir()])


def internal_ip():
    try:
        return socket.gethostbyname(socket.gethostname())
    except socket.gaierror:
        return '127.0.0.1'  # this usually happens on dev laptops; cloud machines work fine


def start_create(sra_id, location):
    if not location.startswith(internal_ip()):
        logger.info(f'Copying from {location}')
        return_code = subprocess.call(['scp', '-r', location, download_dir()])
        if return_code != 0:
            logger.warning(f'Copying from {location} failed')
            return False
        else:
            logger.info(f'Copying from {location} completed successfully')
    else:
        logger.info('Luckily the files to use for building the graph already are on this machine.')
    input_dir = os.path.join(download_dir(), sra_id)
    create_processes[sra_id] = subprocess.Popen(['./create.sh', sra_id, input_dir, create_dir()])
    return True


def start_clean(sra_id, location):
    if not location.startswith(internal_ip()):
        logger.info(f'Copying from {location}')
        return_code = subprocess.call(['scp', '-r', location, args.output_dir])
        if return_code != 0:
            logger.warning('Copying from {location} failed')
            return False
        else:
            logger.info('Copying from {location} completed successfully')
    else:
        logger.info('Luckily the files to use for building the graph already are on this machine.')
    input_file = create_file(sra_id)
    clean_processes[sra_id] = subprocess.Popen(['./clean.sh', sra_id, input_file, clean_dir()])


def is_full():
    return len(download_processes) == MAX_DOWNLOAD_PROCESSES and (
            len(clean_processes) == MAX_CLEAN_PROCESSES or len(create_processes) == MAX_CREATE_PROCESSES)


def ack(operation, sra_id, location):
    url = f'http://{args.server_host}:{args.server_port}/jobs/ack/{operation}'
    data = f'id={sra_id}&location={location}'
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
                logger.warning(f'Server returned response code {response.getcode()} for {url}')
                time.sleep(5)  # avoid overwhelming the server
        except urllib.error.URLError as e:
            print(f'Failed to open URL {url} Reason: {e.reason}')
            time.sleep(5)  # wait a bit and try again


def nack(operation, sra_id):
    url = f'http://{args.server_host}:{args.server_port}/jobs/nack/{operation}'
    data = f'id={sra_id}'
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
                logger.warning(f'Server returned response code {response.getcode()} for {url}')
                time.sleep(5)  # avoid overwhelming the server
        except urllib.error.URLError as e:
            print(f'Failed to open URL {url} Reason: {e.reason}')
            time.sleep(5)  # wait a bit and try again


def check_status():
    completed_downloads = set()
    for sra_id, download_process in download_processes.items():
        return_code = download_process.poll()
        if return_code is not None:
            if return_code == 0:
                logger.info("Download for SRA id {sra_id} completed successfully.")
                completed_downloads.add(sra_id)
                location = os.path.join(f'{internal_ip()}:{download_dir()}', sra_id)

                if not ack('download', sra_id, location):  # all done, yay!
                    return False;
            else:
                logger.warning("Download for SRA id {sra_id} failed.")
                if not nack(sra_id, 'download'):
                    return False  # this shouldn't happen, as we still have work to do
    for d in completed_downloads:
        del download_processes[d]

    completed_creates = set()
    for sra_id, create_process in create_processes.items():
        return_code = create_process.poll()
        if return_code is not None:
            completed_creates.add(sra_id)
            if return_code == 0:
                logger.info("Bulding graph for SRA id {sra_id} completed successfully.")
                location = f'{internal_ip()}:{create_file()}'

                if not ack('create', sra_id, location):  # all done, yay!
                    return False;
            else:
                logger.warning("Building graph for SRA id {sra_id} failed.")
                if not nack(sra_id, 'create'):
                    return False  # this shouldn't happen, as we still have work to do
    for d in completed_creates:
        del create_processes[d]

    completed_cleans = set()
    for sra_id, clean_process in clean_processes.items():
        return_code = clean_process.poll()
        if return_code is not None:
            completed_cleans.add(sra_id)
            if return_code == 0:
                logger.info("Cleaning graph for SRA id {sra_id} completed successfully.")
                location = f'{internal_ip()}:{clean_file()}'

                if not ack('clean', sra_id, location):  # all done, yay!
                    return False;

                # TODO:transfer the cleaned file onto leomed
            else:
                logger.warning("Cleaning graph for SRA id {sra_id} failed.")
                if not nack(sra_id, 'clean'):
                    return False  # this shouldn't happen, as we still have work to do
    for d in completed_cleans:
        del clean_processes[d]

    return True


def do_work():
    while True:
        work_response = get_work()
        print(f'Response is {work_response}')
        if work_response is None:
            break
        if 'download' in work_response:
            start_download(work_response['download']['id'])
        if 'create' in work_response:
            if not start_create(work_response['create']['id'], work_response['create']['location']):
                nack('create', work_response['create']['id'])
        if 'clean' in work_response:
            if not start_clean(work_response['clean']['id'], work_response['clean']['location']):
                nack('clean', work_response['clean']['id'])
        if not check_status():
            break
        time.sleep(5)


def check_env():
    rc = subprocess.call('prereq.sh')
    if rc != 0:
        logger.error("Some pre-requisites are not satisfied, bailing out")
        exit(1)
    pathlib.Path(download_dir()).mkdir(parents=True, exist_ok=True)
    pathlib.Path(create_dir()).mkdir(parents=True, exist_ok=True)
    pathlib.Path(clean_dir()).mkdir(parents=True, exist_ok=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--server_host', help='HTTP server name or ip')
    parser.add_argument('--server_port', default=8000, help='HTTP Port on which the server runs')
    parser.add_argument(
        '--data_dir',
        default=os.path.expanduser('~/Downloads/sra_test/'),
        help='Location of the directory containing the input data')
    parser.add_argument(
        '--output_dir',
        default=os.path.expanduser('~/.metagraph/'),
        help='Location of the directory containing the input data')

    args = parser.parse_args()
    if not os.path.isabs(args.data_dir):
        logger.error(f'data_dir must be an absolute path, not {args.data_dir}')
        exit(1)
    if not args.server_host:
        logger.error('missing --server_host. Can\'t connect to server without it!')
        exit(1)

    do_work()
