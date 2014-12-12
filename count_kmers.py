import pdb
import sys
import gzip
import time
from Bio import SeqIO
import pdb


def print_qgrams():
    f.write('\n\n')
    for key in sorted(hash.keys()):
        f.write(str(key) + ': '+ str(hash[key]) + '\n')




if sys.argv[1] == '-':
    file_path = sys.stdin
elif sys.argv[1].endswith('gz'):
    file_path = gzip.open(sys.argv[1], 'r')
else:
    file_path = open(sys.argv[1], 'r')
q = int(sys.argv[2])
print_qgram_bool = int(sys.argv[3])
myfile = sys.argv[4]

hash = {}
max = 0
mean = 0.0

f = open(myfile, 'w')

t0 = time.time()
for s, seq in enumerate(SeqIO.parse(file_path, "fastq")):
    if s % 100000 == 0:
        t1 = time.time()
        print '%s took %i secs' % (s, t1 - t0)
        t0 = t1

    seq_str = str(seq.seq)
    for i in xrange(0,len(seq_str)-q+1):
        qgram = seq_str[i:i+q]
        try:
            hash[qgram] += 1
        except KeyError:
            hash[qgram] = 1

#print stats

f.write('Number of qgrams: '+ str(len(hash))+'\n')


for key in (hash.keys()):
    mean += hash[key]
    if (hash[key] > max):
        max = hash[key]

mean /= len(hash)

f.write('Maximal appearence of qgram: ' + str(max) + '\nAverage appearence of qgram ' +  str(mean))


if (print_qgram_bool == 1):
    print_qgrams()




f.close()



