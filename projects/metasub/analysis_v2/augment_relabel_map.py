import pandas as pd
import sys
import re

if len(sys.argv) < 3:
    sys.stderr.write('Usage: %s <metadata> <map>\n' % (sys.argv[0]))
    sys.exit(1)

metadata = sys.argv[1] 
infile = sys.argv[2]

### columns of interest
cols = ['metasub_name',
        'city',
        'latitude',
        'longitude',
        'surface_material',
        'station',
        'num_reads',
        'city_latitude',
        'city_longitude',
        'city_total_population',
        'continent',
        'sample_type']

data = pd.read_csv(metadata, delimiter=',', index_col=False)

for line in open(infile, 'r'):
    sl = line.strip(). split('\t')
    extra_fields = ';'.join(['='.join([str(x), str(y)]) for x,y in data[data.uuid == sl[1]][cols].to_dict('records')[0].items()])
    extra_fields = re.sub(' ', '_', extra_fields)
    sl[-1] += ';' + extra_fields
    print('\t'.join(sl))
