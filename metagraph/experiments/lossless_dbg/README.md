# Path Encoder

Building encoder as part of the metagraph library:
```bash
mkdir build
cd build
cmake ../../..
make -j8
make install # to install path_encoder_toolbox
```
When compiled but not installed, executable can be found in `experiments/lossless_dbg/path_encoder_toolbox` and tests in 
`experiments/lossless_dbg/tests`.

Building encoder with global metagraph installation:
```bash
mkdir build
cd build
cmake ..
make -j8
make install # to install path_encoder_toolbox
```

## Examples
```bash
./path_encoder_toolbox compress --input ../tests/data/transcripts_1000.fa --output ./
```

## Query functionality

#### Work in Progress

## Defines
`DEBUG_ASSUME_ALL_KMERS_COVERED` that enables additional assertions under the assumption that for every kmer there is at least one path going through it (covering it)

`DEBUG_ASSUME_NO_TRANSFORMATIONS` enables additional assertions under the assumption that no path transformation is happening

`DEBUG_ADDITIONAL_INFORMATION` more detailed error reports with better insights (beware: slows down the execution as it stores additional information) 

### Unsorted examples

```bash
# logan trim sequences
java -jar trimmomatic-0.39.jar PE -phred33 SRX016044_1.fastq SRX016044_2.fastq  SRX016044_1_filt.fastq SRX016044_1_ufilt.fastq  SRX016044_2_filt.fastq  SRX016044_2_ufilt.fastq ILLUMINACLIP:TruSeq2-PE.fa:2:30:12 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:40 BASECOUNT:N:0:0 

# [SRR554369_1, ERR049156_1, ERR089806_1, SRR327342_1, SRR870667_1, ERR532393_1]
# [SRX016044_1, SRX016044_1_filt, SRX016044_1_ufilt, SRX016044_2, SRX016044_2_filt, SRX016044_2_ufilt]

vsnakemake --profile cluster --dryrun "compressors_results/cwr (kmer-length: 21) (chunks: 1000, statistics-verbosity: 7, path-rerouting: no)/[SRR554369_1, ERR049156_1, ERR089806_1, SRR327342_1, SRR870667_1, ERR532393_1, chr1_10_individuals, chr1_30_individuals]"

vsnakemake --profile benchmark --dryrun "compressors_results/logan/ERR532393_1" 'decompressed_files/cwr/ERR532393_1.fasta' "compressors_results/cwr (kmer-length: 21) (chunks: 0, path-rerouting: no)/ERR532393_1" "compressors_results/cwr (kmer-length: 21) (chunks: 1000, statistics-verbosity: 7, path-rerouting: no)/ERR532393_1"
```

