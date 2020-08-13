import fileinput
import logging
import os
import subprocess
import urllib


def convert_bucket(data_dir):
    """ Converts the original gs://... bucket urls in bucket1..bucket9 to sra_id bucket_no pairs for faster parsing """

    logging.info('Converting bucket files to a faster to parse format...')
    for i in range(1, 10):
        fileout = os.path.join(data_dir, f'bucket_proc{i}')
        with open(fileout, 'w') as fpout:
            file = os.path.join(data_dir, f'bucket{i}')
            try:
                logging.debug(f'Converting {file}')
                with open(file) as fp:
                    for line in fp:
                        parsed = urllib.parse.urlparse(line)
                        sra_id = parsed.path[1:-2]
                        fpout.write(f'{sra_id}\n')
            except FileNotFoundError:
                logging.fatal(f'Could not find {file} in {data_dir}. Please copy it there, otherwise I can\'t '
                              f'figure out which GCloud bucket each sra id is in. Look, I\'m trying to be helpful.')
                exit(1)


class Sra:
    """ Reads a list of files containing SRA ids and returns the ids one by one for further processing """

    def __init__(self, data_files, add_gcloud_bucket, data_dir, ignore_sras):
        self.data_files = data_files
        self.sraid_to_bucket = {}
        self.data_dir = data_dir
        self.ignore_sras = ignore_sras
        if add_gcloud_bucket:
            self.add_buckets()
        self.total_sras = self.get_total_sras()
        self.processed_sras = 0

    def get_total_sras(self):
        logging.info('Computing total number of SRAs...')
        total_sra_count = 0
        ignore_count = 0
        for file in self.data_files:
            with open(file) as fp:
                for line in fp:
                    total_sra_count += 1
                    if line.split(' ')[0] in self.ignore_sras:
                        ignore_count += 1
        to_process_count = total_sra_count - ignore_count
        logging.info(
            f'Found {total_sra_count} SRAs, left out {ignore_count} SRAs, '
            f'processing {to_process_count} SRAs')
        return to_process_count

    def next_item(self):
        """Lazy function (generator) to read a sequence of files line by line."""
        for file in self.data_files:
            with open(file) as fp:
                for line in fp:
                    line = line.rstrip().split()
                    if len(line) == 0:
                        continue
                    if line[0] in self.ignore_sras:  # already processed
                        print(f'{line[0]} is marked as already processed. Skipping...')
                        continue
                    self.processed_sras += 1
                    yield line
        return None

    def load_sraid_to_bucket(self):
        if self.sraid_to_bucket:
            return  # already loaded
        if not all([os.path.exists(os.path.join(self.data_dir, f'bucket{i}')) for i in range(1, 10)]):
            logging.info(
                'SRA Bucket files are not present. Will start downloading. Go get a coffee, this will take a while')
            for i in range(1, 10):
                logging.info(f'Downloading bucket{i}/9...')
                with open(os.path.join(self.data_dir, f'bucket{i}'), 'w') as f:
                    if subprocess.call(['gsutil', '-u', 'metagraph', 'ls', f'gs://sra-pub-run-{i}/'],
                                       stdout=f,
                                       stderr=subprocess.PIPE) != 0:
                        logging.error('Cannot download bucket information from GCS. Heroically dying.')
                        exit(1)
        else:
            logging.info(f'Bucket files found in {self.data_dir}')
        for i in range(1, 10):
            bucket_file = os.path.join(self.data_dir, f'bucket_proc{i}')
            if not os.path.exists(bucket_file):
                convert_bucket(self.data_dir)
            logging.debug(f'Loading {bucket_file} ...')
            with open(bucket_file) as fp:
                for line in fp:
                    self.sraid_to_bucket[line.strip()] = i

    def add_buckets(self):
        """ NCBI data on Google cloud is stored in 9 buckets, named gs://sra-pub-run-1..9
            So for each SRA id, we need to know in which bucket to find it. This function looks at the files in
            args.data_dir and if the ids don't contain a bucket id, it adds it.
        """
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
