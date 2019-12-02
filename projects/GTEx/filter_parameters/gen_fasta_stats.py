import sys
import gzip
import re

lens = []

seq = []
for line in gzip.open(sys.argv[1], 'r'):
    if line[0] == '>':
        if len(seq) > 0:
            _tmp = ''.join(seq)
            lens.append(len(_tmp))
        seq = []
        continue
    seq.append(line.strip())

### get stats
total_len = sum(lens)
total_cnt = len(lens)
median = sorted(lens)[int(total_cnt / 2)]
mean = total_len / float(total_cnt)
max_len = max(lens)

out = open(re.sub(r'.fasta.gz$', '', sys.argv[1]) + '.stats', 'w')
out.write('\t'.join([str(total_cnt), str(total_len), str(median), str(mean), str(max_len)]) + '\n')
out.close()
