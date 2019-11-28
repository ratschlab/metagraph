#!/usr/bin/env python2


import sys
from glob import glob
import os
import commands
import signal
from multiprocessing import Pool


def download(run):
    URL = '/sra/sra-instant/reads/ByRun/sra/{}/{}/{}/{}.sra'.format(run[:3], run[:6], run, run)
    print("Start downloading {}".format(URL))
    output = commands.getstatusoutput("./download_sra_ncbi.sh {}".format(URL))
    print("Finished downloading {}, exit code: {}".format(run, output[0]))
    print(output[1])


def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def main():
    print(sys.argv[0])
    if len(sys.argv) != 3:
        print("Usage: {} <SRA_RUNS_LIST> <num_threads>".format(sys.argv[0]))
        exit(1)

    _, list_filename, num_threads = sys.argv

    runs = []
    with open(list_filename) as f:
        runs = [line.strip() for line in f if len(line.strip())]

    pool = Pool(int(num_threads), init_worker)
    try:
        result = pool.imap_unordered(download, runs)
        pool.close()
        for out in result:
            pass
        print "Finished!"
    except KeyboardInterrupt:
        pool.terminate()
        print "Interrupted!"

    pool.join()


if __name__ == '__main__':
    main()
