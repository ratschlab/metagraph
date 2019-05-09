# Path Encoder

Building encoder as part of the metagraph library:
```bash
mkdir build
cd build
cmake ../../..
make
make install # to install path_encoder_toolbox
```
Executable can be found in `experiments/lossless_dbg/path_encoder_toolbox` and tests in 
`experiments/lossless_dbg/tests`.

Building encoder with global metagraph installation:
```bash
mkdir build
cd build
cmake ..
make
make install # to install path_encoder_toolbox
```

## Examples
```bash
path_encoder_toolbox compress --compressor-type wavelet --statistics statistics.json --input SRR554369_1.fasta --output SRR554369_1/
```