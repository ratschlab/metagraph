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

def _print(curr_r1, curr_r2):
    global out_r1
    global out_r2
    global maxlen
    print('@' + curr_r1[0] + '/1', file=out_r1)
    print(curr_r1[2][:maxlen], file=out_r1)
    print('+', file=out_r1)
    print(curr_r1[1][:maxlen], file=out_r1)
    print('@' + curr_r2[0] + '/2', file=out_r2)
    print(curr_r2[2][:maxlen], file=out_r2)
    print('+', file=out_r2)
    print(curr_r2[1][:maxlen], file=out_r2)


samfile = pysam.AlignmentFile(fname, 'rb')

out_r1 = open(re.sub(r'.bam$', '', fname) + '.r1.fq', 'w') 
out_r2 = open(re.sub(r'.bam$', '', fname) + '.r2.fq', 'w')

last_id = ''
curr_r1 = []
curr_r2 = []
cnt = 0
warn_count = 0

for read in samfile.fetch(until_eof=True):
    if cnt > 0 and cnt % 100000 == 0:
        print('[processed %.2f million read pairs]' % (cnt / 1000000.0), file=sys.stderr)
    ### ignore secondary
    if read.is_secondary:
        continue
    curr_id = read.qname
    #if len(last_id) > 0 and last_id > curr_id:
    #    print >> sys.stderr, 'Order of read IDs locally decreasing - BAM file seems not to be sorted by read ID'
    #    sys.exit(1)
    if last_id != curr_id:
        if len(last_id) > 0:
            if len(curr_r1) == 0:
                if warn_count < 10:
                    print('WARNING: ' + curr_r2[0] + ' found to be unpaired', file=sys.stderr)
                warn_count += 1
            elif len(curr_r2) == 0:
                if warn_count < 10:
                    print('WARNING: ' + curr_r1[0] + ' found to be unpaired', file=sys.stderr)
                warn_count += 1
            else:
                _print(curr_r1, curr_r2)
                cnt += 1
        last_id = curr_id

    if read.is_reverse:
        curr_r = [read.qname, read.qual[::-1], _rev(read.seq)]
    else:
        curr_r = [read.qname, read.qual, read.seq]
    if read.is_read1:
        curr_r1 = curr_r 
    else:
        curr_r2 = curr_r

if len(curr_r1) == 0:
    if warn_count < 10:
        print(sys.stderr, 'WARNING: ' + curr_r2[0] + ' found to be unpaired')
    warn_count += 1
elif len(curr_r2) == 0:
    if warn_count < 10:
        print(sys.stderr, 'WARNING: ' + curr_r1[0] + ' found to be unpaired')
    warn_count += 1
else:
    _print(curr_r1, curr_r2)

out_r1.close()
out_r2.close()
            
if warn_count >= 10:
    print(sys.stderr, 'SUPPRESSED %i Warnings' % (warn_count - 10))
