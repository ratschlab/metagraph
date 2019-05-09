# Path Encoder

Building encoder as part of the metagraph library:
```bash
mkdir build
cd build
cmake ../../..
make
make install # to install path_encoder_toolbox
```
When not installed, executable can be found in `experiments/lossless_dbg/path_encoder_toolbox` and tests in 
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

path_encoder_toolbox compress --compressor-type wavelet --statistics statistics.json --input INPUT.fasta --output OUTPUT_FOLDER/
```