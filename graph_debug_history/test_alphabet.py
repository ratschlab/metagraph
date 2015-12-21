import sys

if len(sys.argv) < 2:
    print >> sys.stderr, "Usage: %s <fasta.fa>" % sys.argv[0]
    sys.exit(1)

alph = dict()
for l, line in enumerate(open(sys.argv[1], 'r')):
    if l > 0 and l % 100000 == 0:
        sys.stdout.write('.')
        if l % 1000000 == 0:
            sys.stdout.write('%i' % l)
        sys.stdout.flush()
    if line[0] == '>':
        continue
    for c in line.strip():
        try:
            alph[c] += 1
        except KeyError:
            alph[c] = 1

for k, v in alph.iteritems():
    print "%s : %i" % (k, v)
