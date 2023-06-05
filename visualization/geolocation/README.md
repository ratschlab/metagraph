# Search Geo Coordinates for DNA sequences

#### This front end part was moved to a separate repository https://github.com/ratschlab/projects2019-dna_web_search

1. Compile metagraph
2. Run this
```
exe=../../metagraph/build/metagraph

input_fasta="../../metagraph/tests/data/transcripts_1000.fa"
geo_fasta="data/random_seq_locations.fa"

python3 data/generate_random_annotation.py $input_fasta > $geo_fasta

$exe build -k 10 -o data/test_graph $geo_fasta
$exe annotate --anno-header -i data/test_graph -o data/test_annotation $geo_fasta

python2 sequence_annotator.py data/test_graph data/test_annotation.column.annodbg > data/random_locations.js

open world_data.html
```

## Test Server

1. Run server
```
../../metagraph/build/metagraph server_query -i data/test_graph -a data/test_annotation.column.annodbg
```
2. Run client
```
python3 client_annotator.py localhost 5555 > data/random_locations.js
```
3. Check results
```
open world_data.html
```
