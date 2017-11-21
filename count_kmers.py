import pdb
import sys
import gzip
import time
from Bio import SeqIO
import pdb
import itertools
import collections
import numpy as np


if sys.argv[1] == '-':
    file_path = sys.stdin
elif sys.argv[1].endswith('gz'):
    file_path = gzip.open(sys.argv[1], 'r')
else:
    file_path = open(sys.argv[1], 'r')
q = int(sys.argv[2])
print_qgram = bool(sys.argv[3])
myfile = sys.argv[4]

counter = collections.Counter()

f = open(myfile, 'w')

t0 = time.time()
for fmt in "fastq", "fasta":
    file_path.seek(0)
    for s, seq in enumerate(SeqIO.parse(file_path, fmt)):
        if s % 100000 == 0:
            t1 = time.time()
            print('%s took %i secs' % (s, t1 - t0))
            t0 = t1

        seq_str = str(seq.seq)
        counter.update((seq_str[i:i+q] for i in range(len(seq_str) - q + 1)))

f.write('Number of qgrams: ' + str(len(counter)) + '\n'
        'Maximal appearence of qgram: ' + str(max(counter.values())) + '\n'
        'Average appearence of qgram ' +  str(np.mean(list(counter.values()))))

if print_qgram:
    f.write('\n\n')
    for key in sorted(counter.keys()):
        f.write(str(key) + ': '+ str(counter[key]) + '\n')

f.close()
