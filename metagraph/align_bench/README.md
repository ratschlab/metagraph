# Default vs. redone aligner comparison

Benchmark-only tooling (branch `bench/align-engine-compare`). Nothing here is
built by default.

- `metagraph align --align-engine redone` runs the opt-in engine
  (`graph/alignment_redone`) from the CLI: TSV only, best path only.
- `align_engine_bench` (CMake target, `EXCLUDE_FROM_ALL`) runs one engine per
  process on the same graph and annotation. It times each read, has a per-read
  watchdog, and also has a `fanout` mode. Usage is in the source header.

## Recipe

Run everything from a work directory `W` that holds `data/` and `runs/`.

1. Build: `make align_engine_bench metagraph` in a Release build directory.
2. Genomes: 24 RefSeq "reference genome" / "Complete Genome" bacteria. The
   accessions are in `genomes.tsv` next to this file, and each genome is one
   FASTA file `data/genomes/<accession>.fa`. Download them from the
   `ftp_path` column of RefSeq's `bacteria/assembly_summary.txt`.
3. Graphs, k=31 (`M` = the metagraph binary):
   ```
   M build -k 31 --mode canonical -p 16 -o data/graphs/canonical data/genomes/*.fa
   M build -k 31 --mode basic     -p 16 -o data/graphs/basic     data/genomes/*.fa
   M transform --to-fasta --primary-kmers -p 16 -o data/graphs/primary_contigs data/graphs/canonical.dbg
   M build -k 31 --mode primary -p 16 -o data/graphs/primary data/graphs/primary_contigs.fasta.gz
   for m in canonical basic primary; do
     M annotate -i data/graphs/$m.dbg --anno-filename --coordinates -p 16 -o data/graphs/$m data/genomes/*.fa
   done
   ```
   This gives one label per genome file and a column annotation with
   coordinates.
4. Reads (badread from a venv; seeds are fixed):
   - `data/reads/allref.fa` = `cat data/genomes/*.fa`.
   - `data/reads/contig_map.tsv` has the columns contig, genome and contig
     length.
   - Then run `simulate.sh` (next to this file), from `data/`:
     - Illumina-like reads: 2000 × 150 bp, 0.5% substitutions, + 100 random
       negatives (`simulate_illumina.py`).
     - ONT-like reads (badread `nanopore2023`, no adapters, junk, random
       reads, chimeras or glitches) at identity 97/95/90%:
       - 1000 × 300–500 bp (mean 400);
       - 300 × 1–5 kb (mean ~2.9 kb).
     - Each ONT set gets +5% random negatives (`badread_to_fasta.py`).
   - Truth is in each read name: `r<i>|<genome>|<contig>|<forward start>|<strand>`.
5. Matrix: `run_all.sh <bench binary> 4`. Each read gets a 10 s cap, and an
   arm stops after 25 timed-out reads.
6. Scoring:
   - `metrics.py data/genomes data/reads/<set>.fa runs/<arm>/<engine>.<set> ...`
   - `fanout_summary.py runs/fanout/*.tsv`, after
     `align_engine_bench fanout 0 0 align -i data/graphs/canonical.dbg --align-min-seed-length 19 <reads>`.
