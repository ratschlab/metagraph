import argparse
import json
import logging
import os
import time
import urllib.request

logger = logging.getLogger('metagraph-client')
args = None
download_count = 0
create_count = 0
clean_count = 0


def get_work():
    url = f'http://{args.server_host}:{args.server_port}/jobs?download={download_count}&create={create_count}&clean={clean_count}'
    for i in range(10):
        try:
            response = urllib.request.urlopen(url)
            if response.getcode() == 200:
                return json.loads(response.read())
            else:
                logger.warning(f'Server returned response code {response.getcode()} for {url}')
                time.sleep(5)  # avoid overwhelming the server
        except urllib.error.URLError as e:
            print(f'Failed to open URL {url} Reason: {e.reason}')
            time.sleep(5)  # wait a bit and try again


def do_work():
    is_done = False
    while not is_done:
        work_response = get_work()
        print(f'Response is {work_response}')
        time.sleep(10)


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

    if not args.server_host:
        logger.error('missing --server_host. Can\'t connect to server without it!')
        exit(1)

    do_work()
