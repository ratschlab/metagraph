# given a list of full paths of SRAs stored in GCS, keep the ones that have an id present in #to_copy_name
source_name = '/home/ddanciu/gcs_fungi_files'
to_copy_name = '/home/ddanciu/sra_fungi_metadata_no_pacbio_nanopore_with_threshold.ids'
dest = open('/tmp/to_copy', 'w')

to_copy = set()
lines = open(to_copy_name).readlines()
for line in lines:
    to_copy.add(line[:-1])
print(f'To copy has {len(to_copy)} elements')
for line in open(source_name):
    sra_id = line.split('/')[-1].split('.')[0]
    if sra_id in to_copy:
        to_copy.remove(sra_id)
        dest.write(line)

