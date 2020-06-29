#!/usr/bin/env python3

import argparse
import collections
import http
import http.server
import json
import logging
import math
import multiprocessing
import os
import psutil
import signal
import subprocess
import threading
import time
import util
import urllib.parse
import urllib.request

args = None

sra_info = {}  # global information about a processed SRA (for now, time only)

pending_processes = []  # TODO: probably not needed
download_processes = {}
build_processes = {}
clean_processes = {}
transfer_processes = {}

waiting_builds = collections.OrderedDict()
waiting_cleans = collections.OrderedDict()

CORES = multiprocessing.cpu_count()
MAX_DOWNLOAD_PROCESSES = CORES / 4
MAX_BUILD_PROCESSES = CORES / 4
MAX_CLEAN_PROCESSES = CORES / 4

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
    if downloads_done or len(download_processes) >= MAX_DOWNLOAD_PROCESSES or len(
            waiting_builds) + len(download_processes) >= 2 * MAX_BUILD_PROCESSES:
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


def destination_dir(sra_id, subdir):
    result = ''
    if len(sra_id) > 6:
        result = f'{sra_id[0:6]}/{result}'
    if len(sra_id) > 3:
        result = f'{sra_id[0:3]}/{result}'
    return os.path.join(args.destination, f'{subdir}/{result}')


def dump_pending(sra_ids):
    ids = ','.join(sra_ids)
    url = f'http://{args.server}/jobs/preempt'
    data = f'client_id={util.internal_ip()}&ids={ids}'
    with open('/tmp/shutdown.sh', 'w') as f:
        f.write(f'curl --data "{data}&reason=shutdown" "{url}"')


def start_download(download_resp):
    if 'id' not in download_resp:
        logging.info('No more downloads available. We\'re almost done!')
        return
    sra_id = download_resp['id']
    pending_processes.append(sra_id)
    dump_pending(pending_processes)
    util.make_dir_if_needed(download_dir(sra_id))
    log_file_name = os.path.join(download_dir(sra_id), 'download.log')
    log_file = util.TeeLogger(log_file_name, 'Stage')
    bucket = '0'
    if args.source == 'ncbi':
        if 'bucket' not in download_resp:
            logging.info(f'[{sra_id}] Specified NCBI as download source, but server response has no "bucket" field. '
                         f'Will download via HTTP instead of GCS')
        else:
            bucket = download_resp['bucket']
    download_processes[sra_id] = (
        subprocess.Popen(['./download.sh', args.source, sra_id, download_dir_base(), bucket], stdout=log_file,
                         stderr=log_file), time.time())
    sra_info[sra_id] = (time.time(),)


def write_log_header(log_file, operation, sra_id, required_ram_gb, available_ram_gb):
    free_ram_gb = psutil.virtual_memory().available / 1e9
    log_file.write(
        f'[{sra_id}] Starting {operation} on {util.internal_ip()}, required RAM {round(required_ram_gb, 2)}, free RAM '
        f'{round(free_ram_gb, 2)}GB, available for {operation} (not reserved) RAM {round(available_ram_gb, 2)}GB')
    log_file.write(f'[{sra_id}] Full machine log: "gsutil cat {args.destination}logs/{util.internal_ip()}/client.log"')


def start_build(sra_id, wait_time, buffer_size_gb, container_type, required_ram_gb, available_ram_gb):
    input_dir = download_dir(sra_id)
    output_dir = build_dir(sra_id)
    num_cores = max(4, round(required_ram_gb / 3.75))
    logging.info(f'[{sra_id}] Starting build from {input_dir} to {output_dir}, buffer={round(buffer_size_gb, 2)}GB '
                 f'on {num_cores} cores')
    util.make_dir_if_needed(build_dir(sra_id))
    log_file_name = os.path.join(build_dir(sra_id), 'build.log')
    log_file = util.TeeLogger(log_file_name)
    write_log_header(log_file, 'build', sra_id, required_ram_gb, available_ram_gb)
    build_processes[sra_id] = (subprocess.Popen(
        ['./build.sh', sra_id, input_dir, output_dir, str(buffer_size_gb), str(num_cores), container_type], stdout=log_file,
        stderr=log_file), time.time(), wait_time, required_ram_gb)
    return True


def start_clean(sra_id, wait_time, kmer_count_singletons, fallback, required_ram_gb, available_ram_gb):
    input_file = build_file(sra_id)
    output_file = os.path.join(clean_dir(sra_id), sra_id)
    num_cores = max(4, round(required_ram_gb / 3.75))
    logging.info(f'[{sra_id}] Starting clean from {input_file} to {output_file} on {num_cores} cores')
    log_file_name = os.path.join(clean_dir(sra_id), 'clean.log')
    log_file = util.TeeLogger(log_file_name, '%,')  # %, eliminates the progress bar spam
    write_log_header(log_file, 'clean', sra_id, required_ram_gb, available_ram_gb)
    clean_processes[sra_id] = (
        subprocess.Popen(
            ['./clean.sh', sra_id, input_file, output_file, str(kmer_count_singletons), str(fallback), str(num_cores)],
            stdout=log_file,
            stderr=log_file),
        time.time(), wait_time, required_ram_gb)
    return True


def start_transfer(sra_id, source, top_folder):
    destination = destination_dir(sra_id, top_folder)
    transfer_process = subprocess.Popen(
        f'gsutil -q -u metagraph cp -r {source} {destination}', shell=True), time.time()
    if top_folder == 'clean':
        transfer_processes[sra_id] = transfer_process
    if sra_id in pending_processes:
        pending_processes.remove(sra_id)
        dump_pending(pending_processes)
    else:
        logging.error(f'[sra_id] just transferred successfully but not present in pending_ids. Something is messed up.')


def ack(operation, params):
    url = f'http://{args.server}/jobs/ack/{operation}'
    data = f'client_id={util.internal_ip()}&' + urllib.parse.urlencode(params)
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
    data = f'client_id={util.internal_ip()}&' + urllib.parse.urlencode(params)
    sra_id = params['id']
    cleaned_dir = clean_dir(sra_id)
    start_transfer(sra_id, cleaned_dir, 'not_' + operation)
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


# Simple sanity check for graph built: unique kmer count from KMC * 2 is equal to edges - dummy_source - dummy_sink
def check_sanity(sra_id):
    stat_file = open(os.path.join(build_dir(sra_id), f'{sra_id}.stats')).readlines()
    if len(stat_file) != 3:
        logging.error(f'[{sra_id}] Stat file should have 3 lines, but has ${len(stat_file)}')
        return False
    unique_kmers = int(stat_file[0].split(' ')[-1])
    unique_kmers -= (int(stat_file[1].split(' ')[-1]) + int(stat_file[2].split(' ')[-1]))
    sanity = unique_kmers == 2 * sra_info[sra_id][2]  # *2 because we build a canonical graph
    if not sanity:
        logging.error(
            f'[{sra_id}] Sanity check for graph build failed. Expected 2 * {sra_info[sra_id][2]} kmers, got {unique_kmers}')
    return sanity


def check_status():
    global must_quit
    if must_quit:
        return False
    completed_downloads = set()
    for sra_id, (download_process, start_time) in download_processes.items():
        return_code = download_process.poll()
        is_timed_out = (time.time() - start_time) > 120 * 60
        if return_code is not None or is_timed_out:
            if os.path.exists(os.path.join(download_dir(sra_id), 'code')):
                return_code = int(open(os.path.join(download_dir(sra_id), 'code')).read())
            elif is_timed_out:
                logging.warning(f'[{sra_id}] Download timed out after {time.time() - start_time} seconds.')
                return_code = 254
            else:
                logging.error(f'[{sra_id}] Download process did not provide a return code. Assuming error')
                return_code = 255
            completed_downloads.add(sra_id)
            log_file_name = os.path.join(download_dir(sra_id), 'download.log')
            logging.info(f'[{sra_id}] Download finished with output\n {util.tab(open(log_file_name).readlines())}\n\n\n')
            log_new_file_name = os.path.join(clean_dir(sra_id), 'download.log')
            util.make_dir_if_needed(clean_dir(sra_id))
            os.rename(log_file_name, log_new_file_name)

            download_path = download_dir(sra_id)
            sra_dir = os.path.join(download_path, 'sra')
            download_size_mb = util.dir_size_MB(sra_dir)
            size_file = os.path.join(download_dir(sra_id), 'size')
            if os.path.exists(size_file):
                sra_size_mb = round(int(open(size_file).read()) / 1e6, 2)
                logging.info(f'Downloaded sra files have {sra_size_mb}MB')
            else:
                logging.warning('Could not find size file. Reporting size -1')
                sra_size_mb = -1
            subprocess.run(['rm', '-rf', sra_dir])
            kmc_dir = os.path.join(download_path, 'kmc')
            kmc_size_mb = util.dir_size_MB(kmc_dir)
            success = True
            if return_code == 0:
                logging.info(f'[{sra_id}] Download completed successfully.')
                stats_file = os.path.join(download_path, 'stats')
                try:
                    with open(stats_file) as stats:
                        json_resp = json.loads(stats.read())
                    if '#k-mers_coverage' in json_resp and '#k-mers_below_min_threshold' in json_resp:
                        kmer_count_unique = int(json_resp['#Unique_counted_k-mers'])
                        kmer_coverage = int(json_resp['#k-mers_coverage'])
                        kmer_count_singletons = int(json_resp['#k-mers_below_min_threshold'])
                    else:
                        logging.warning(f'[{sra_id}] Invalid KMC stat files, assuming failure')
                        success = False
                except FileNotFoundError:
                    logging.warning(f'[{sra_id}] Could not find KMC stats file {stats_file}, baling out.')
                    success = False
            else:
                success = False
            if success:
                params = {'id': sra_id, 'time': int(time.time() - start_time), 'size_mb': sra_size_mb,
                          'download_size_mb': download_size_mb, 'kmc_size_mb': kmc_size_mb,
                          'kmer_count_unique': kmer_count_unique, 'kmer_coverage': kmer_coverage,
                          'kmer_count_singletons': kmer_count_singletons}
                sra_info[sra_id] = (
                    *sra_info[sra_id], sra_size_mb, kmer_count_unique, kmer_coverage, kmer_count_singletons)
                ack('download', params)
                waiting_builds[sra_id] = (time.time())
            else:
                logging.warning(f'[{sra_id}] Download failed. Removing {download_path}')
                subprocess.run(['rm', '-rf', download_path])
                params = {'id': sra_id, 'time': int(time.time() - start_time), 'size_mb': sra_size_mb,
                          'download_size_mb': download_size_mb, 'kmc_size_mb': kmc_size_mb, 'exit_code': return_code}
                nack('download', params)
    for d in completed_downloads:
        del download_processes[d]

    completed_builds = set()
    total_reserved_ram_gb = 0  # how much memory all active processes need
    used_cores = 0
    for sra_id, (build_process, start_time, wait_time, reserved_ram_gb) in build_processes.items():
        return_code = build_process.poll()

        if return_code is not None:
            completed_builds.add(sra_id)
            log_file_name = os.path.join(build_dir(sra_id), 'build.log')
            logging.info(f'[{sra_id}] Build finished with output\n {util.tab(open(log_file_name).readlines())}\n\n\n')
            log_new_file_name = os.path.join(clean_dir(sra_id), 'build.log')
            os.rename(log_file_name, log_new_file_name)

            # clean up the download path; if adding retries, do this only on success
            download_path = download_dir(sra_id)
            logging.info(f'[{sra_id}] Cleaning up {download_path}')
            subprocess.run(['rm', '-rf', download_path])

            build_path = build_dir(sra_id)
            build_size_mb = util.dir_size_MB(build_path)
            if return_code == 0:
                logging.info(f'[{sra_id}] Building graph completed successfully.')
                sanity = check_sanity(sra_id)
                params = {'id': sra_id, 'time': int(time.time() - start_time),
                          'wait_time': int(wait_time), 'size_mb': build_size_mb, 'sanity': sanity}
                ack('build', params)
                waiting_cleans[sra_id] = (time.time())
            else:
                logging.warning(f'[{sra_id}] Building graph failed. Removing {build_path}.')
                subprocess.run(['rm', '-rf', build_path])
                params = {'id': sra_id, 'time': int(time.time() - start_time), 'wait_time': int(wait_time),
                          'size_mb': build_size_mb, 'return_code': return_code}
                nack('build', params)
        else:
            total_reserved_ram_gb += reserved_ram_gb
            used_cores += max(4, round(reserved_ram_gb / 3.75))
    for d in completed_builds:
        del build_processes[d]

    completed_cleans = set()
    for sra_id, (clean_process, start_time, wait_time, reserved_ram_gb) in clean_processes.items():
        return_code = clean_process.poll()
        if return_code is not None:
            completed_cleans.add(sra_id)
            log_file_name = os.path.join(clean_dir(sra_id), 'clean.log')
            logging.info(f'[{sra_id}] Clean finished with output\n {util.tab(open(log_file_name).readlines())}\n\n\n')

            # clean up the original graph; if adding retries, do this only on success
            build_path = build_dir(sra_id)
            logging.info(f'[{sra_id}] Cleaning up {build_path}')
            subprocess.run(['rm', '-rf', build_path])

            cleaned_dir = clean_dir(sra_id)
            cleaned_size_mb = util.dir_size_MB(cleaned_dir)
            if return_code == 0:
                logging.info(f'[{sra_id}] Cleaning graph completed successfully.')

                params = {'id': sra_id, 'time': int(time.time() - start_time),
                          'size_mb': cleaned_size_mb}
                ack('clean', params)
                start_transfer(sra_id, cleaned_dir, 'clean')
            else:
                logging.warning(f'[{sra_id}] Cleaning graph failed. Removing {cleaned_dir}')
                subprocess.run(['rm', '-rf', cleaned_dir])
                params = {'id': sra_id, 'time': int(time.time() - start_time),
                          'size_mb': cleaned_size_mb, 'return_code': return_code}
                nack('clean', params)
        else:
            total_reserved_ram_gb += reserved_ram_gb
            used_cores += max(4, round(reserved_ram_gb / 3.75))
    for d in completed_cleans:
        del clean_processes[d]

    completed_transfers = set()
    for sra_id, (transfer_process, start_time) in transfer_processes.items():
        return_code = transfer_process.poll()
        if return_code is not None:
            completed_transfers.add(sra_id)
            # clean up the cleaned graph; if adding retries, do this only on success
            clean_path = clean_dir(sra_id)
            cleaned_size_mb = util.dir_size_MB(clean_path)
            logging.info(f'[{sra_id}] Cleaning up {clean_path}')
            subprocess.run(['rm', '-rf', clean_path])

            if return_code == 0:
                logging.info(f'[{sra_id}] Transferring graph completed successfully.')
                params = {'id': sra_id, 'time': int(time.time() - start_time),
                          'total_time': int(time.time() - sra_info[sra_id][0]), 'size_init_mb': sra_info[sra_id][1],
                          'size_final_mb': cleaned_size_mb}
                ack('transfer', params)
            else:
                logging.warning(f'[{sra_id}] Transferring cleaned graph failed.')
                params = {'id': sra_id, 'time': int(time.time() - start_time),
                          'size_mb': cleaned_size_mb}
                nack('transfer', params)

    # for cleaning we allow using all the available RAM
    total_ram_gb = psutil.virtual_memory().total / 1e9
    available_ram_gb = total_ram_gb - total_reserved_ram_gb
    if used_cores < CORES and waiting_cleans:
        logging.info(f'Ram reserved {round(total_reserved_ram_gb, 2)}GB, total {round(total_ram_gb, 2)}')
        for sra_id, (start_time) in waiting_cleans.items():
            # remove the old clean waiting and append the new one after
            build_path = build_dir(sra_id)
            build_size_gb = util.dir_size_MB(build_path) / 1e3
            required_ram_gb = max(build_size_gb * 1.1, build_size_gb + 1)
            if available_ram_gb > required_ram_gb:
                logging.info(
                    f'[{sra_id}] Estimated {required_ram_gb}GB needed for cleaning, available {available_ram_gb} GB')
                kmer_count_unique = sra_info[sra_id][2]
                kmer_coverage = sra_info[sra_id][3]
                kmer_count_singletons = sra_info[sra_id][4]
                fallback = 5 if kmer_coverage > 5 else 2 if kmer_coverage > 2 or kmer_count_unique > 1e6 else 1

                # multiplying singletons by 2 bc we compute canonical graph and KMC doesn't
                start_clean(sra_id, time.time() - start_time, 2 * kmer_count_singletons, fallback, required_ram_gb,
                            available_ram_gb)
                available_ram_gb -= required_ram_gb
                del waiting_cleans[sra_id]
                break
            logging.info(f'[{sra_id}] Not enough RAM for cleaning. '
                         f'Have {round(available_ram_gb, 2)}GB need {round(build_size_gb + 0.5, 2)}GB')

    if used_cores < CORES and waiting_builds:
        logging.info(f'Ram reserved {round(total_reserved_ram_gb, 2)}GB, total {round(total_ram_gb, 2)}')
        for sra_id, (start_time) in waiting_builds.items():
            num_kmers = sra_info[sra_id][2]
            # estimate RAM needed for loading graph in memory;
            bytes_per_kmer = 2.6  # 0.6 bytes/kmer (for --small representation), 2 byte/kmer-count
            kmer_count = 2.6 * num_kmers  # 2x canonical+non-canonical +  ~30% for dummy kmers (typically it's 10%)
            required_ram_gb = round(kmer_count * bytes_per_kmer / 1e9 + 0.5, 2)
            if required_ram_gb > total_ram_gb - 2:
                download_path = download_dir(sra_id)
                logging.warning(
                    f'[{sra_id}] Building graph needs too much RAM: {required_ram_gb}GB). Removing {download_path}.')
                subprocess.run(['rm', '-rf', download_path])
                params = {'id': sra_id, 'time': int(time.time() - start_time), 'required_ram_gb': required_ram_gb}
                nack('build', params)
                del waiting_builds[sra_id]
                break
            elif required_ram_gb < available_ram_gb and available_ram_gb > 2:
                logging.info(
                    f'[{sra_id}] Estimated {required_ram_gb}GB needed for building, available {available_ram_gb} GB')
                # how much memory does it take to load all unique kmers into RAM: 8B for the kmer, 2B for the count
                required_ram_all_mem_gb = num_kmers * (8 + 2) * 3.5 / 1e9;  # also account for dummy kmers
                if required_ram_all_mem_gb < 5 and required_ram_all_mem_gb < available_ram_gb:
                    required_ram_gb = max(required_ram_gb, required_ram_all_mem_gb)
                    start_build(sra_id, time.time() - start_time, math.ceil(required_ram_all_mem_gb), 'vector',
                                required_ram_gb, available_ram_gb)
                else:
                    buffer_size_gb = max(2, min(round(required_ram_gb * 0.8 - 1), 20))
                    start_build(sra_id, time.time() - start_time, buffer_size_gb, 'vector_disk', required_ram_gb,
                                available_ram_gb)
                del waiting_builds[sra_id]
                available_ram_gb -= required_ram_gb  # not that it matters
                break
            else:
                logging.info(
                    f'[{sra_id}] Not enough RAM for building. Have {round(total_ram_gb - total_reserved_ram_gb, 2)}GB '
                    f'need {required_ram_gb}GB')

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
            logging.error(f'[{sra_id}] Server response invalid. Expected a \'download\' tag: {work_response}')
        subprocess.Popen(
            [f'gsutil rsync -x \'(?!^client.log$)\' {args.output_dir} {args.log_destination}/{util.internal_ip()}/'],
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)

        i = i + 1
        time.sleep(5)
    subprocess.Popen(
        [f'gsutil rsync -x \'(?!^client.log$)\' {args.output_dir} {args.log_destination}/{util.internal_ip()}/'],
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)


def check_env():
    """ Make sure all the necessary software is in place to successfully run the client and create working
    directories """

    util.make_dir_if_needed(download_dir_base())
    util.make_dir_if_needed(build_dir_base())
    util.make_dir_if_needed(clean_dir_base())

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
    data = f'client_id={util.internal_ip()}&ids={ids}'
    try:
        request = urllib.request.Request(url, data=data.encode('UTF-8'),
                                         headers={'Content-type': 'application/x-www-form-urlencoded'},
                                         method='POST')
        response = urllib.request.urlopen(request, timeout=10)
        if response.getcode() != 200:
            logging.warning(f'Server returned response code {response.getcode()} for {url}')
            time.sleep(5)  # avoid overwhelming the server
    except urllib.error.URLError as e:
        logging.error(f'Failed to open URL {url} Reason: {e.reason}')
    except http.client.RemoteDisconnected as e:
        logging.error(f'Failed to open URL {url} Reason: remote disconnected.')
    except:
        logging.error(f'Failed to open URL {url} Reason: unknown.')
    global pending_processes
    pending_processes = []
    dump_pending(pending_processes)
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
    parser.add_argument('--destination', default=None,
                        help='Host/directory where the cleaned BOSS graphs are copied to')
    parser.add_argument('--log_destination', default=None,
                        help='GS folder where client logs are collected')
    parser.add_argument('--port', default=8001, help='HTTP Port on which the status/kill server runs')
    parser.add_argument('--server_info', default=None,
                        help='Where to obtain the server host/port on gcs')
    args = parser.parse_args()

    args.log_destination = args.log_destination or os.path.join(args.destination, 'logs')
    args.server_info = args.server_info or os.path.join(args.destination, 'server')

    if not os.path.isabs(args.output_dir):
        logging.error(f'output_dir must be an absolute path, not {args.output_dir}')
        exit(1)

    check_env()

    logging.info(f'Machine has {CORES} cores. Max simultaneous downloads: {MAX_DOWNLOAD_PROCESSES}. Max builds: '
                 f'{MAX_BUILD_PROCESSES}. Max cleans: {MAX_CLEAN_PROCESSES}')

    if not args.server:
        logging.info('Trying to find server address...')
        if subprocess.call(['gsutil', 'cp', args.server_info, '/tmp/server'], stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE) != 0:
            logging.error('Cannot find server ip/port on Google Cloud Storage. Sorry, I tried.')
            exit(1)
        with open('/tmp/server') as fp:
            args.server = fp.read()

        logging.info(f'Found server at {args.server}')

    # gracefully handle termination
    signal.signal(signal.SIGTERM, signal_handler)
    # signal.signal(signal.SIGKILL, signal_handler)
    signal.signal(signal.SIGSEGV, signal_handler)
    signal.signal(signal.SIGINT, signal_handler)

    httpd = http.server.HTTPServer(('', args.port), SimpleHTTPRequestHandler)
    logging.info(f'Starting client status server on port {args.port}')
    thread = threading.Thread(name='server_thread', target=httpd.serve_forever)
    thread.daemon = True
    thread.start()

    do_work()
