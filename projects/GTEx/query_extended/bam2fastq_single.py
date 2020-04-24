import sys
import pysam
import re

if len(sys.argv) < 2:
    print('Usage: %s <bam_file> [<max_len>]' % sys.argv[0], file=sys.stderr)
    sys.exit(1)
fname = sys.argv[1]
maxlen = 100000
if len(sys.argv) > 2:
    maxlen = int(sys.argv[2])

rc_dict = {'A':'T',
           'T':'A',
           'C':'G',
           'G':'C',
           'N':'N'}

def _rev(s):
    global rc_dict
    return ''.join(rc_dict[_] for _ in s[::-1])

def _print(curr_r):
    global out_r1
    global maxlen
    print('@' + curr_r[0] + '/1', file=out_r1)
    print(curr_r[2][:maxlen], file=out_r1)
    print('+', file=out_r1)
    print(curr_r[1][:maxlen], file=out_r1)

samfile = pysam.AlignmentFile(fname, 'rb')

out_r1 = open(re.sub(r'.bam$', '', fname) + '.r1.fq', 'w') 

last_id = ''
curr_r = []
cnt = 0

for read in samfile.fetch(until_eof=True):
    if cnt > 0 and cnt % 100000 == 0:
        print('[processed %.2f million reads]' % (cnt / 1000000.0), file=sys.stderr)
    curr_id = read.qname
    if last_id != curr_id:
        if len(last_id) > 0:
            _print(curr_r)
            cnt += 1
        last_id = curr_id

    if read.is_reverse:
        curr_r = [read.qname, read.qual[::-1], _rev(read.seq)]
    else:
        curr_r = [read.qname, read.qual, read.seq]

_print(curr_r)

out_r1.close()
