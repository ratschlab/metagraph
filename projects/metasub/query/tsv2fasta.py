import sys

infname = sys.argv[1]

def wprint(s, w=80):
    for i in range(0, len(s), w):
        print s[i:min(i+w, len(s))]

for line in open(infname, 'r'):
    sl = line.strip().split('\t')
    if sl[6] != '*':
        print('>%s' % sl[0])
        wprint(sl[3])
