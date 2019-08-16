#/usr/bin/env bash

DATA_PATH="/cluster/work/grlab/projects/metagenome/data/small_kingsford"
DATA="$DATA_PATH/*[^s].fasta.gz"

NUM_THREADS=5

if [ -f metagraph ]; then
  exe="./metagraph"
else
  exe="$(dirname ${BASH_SOURCE[0]})/../build/metagraph"
fi

echo Test $exe
echo ""

big_file="$(ls -lhS $DATA | head -n 1 | awk '{print $9}')"

time="gtime -f Time...%E\nMem....%M.kB\nCPU....%P\n"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "~~~~~~~~~~~~~~~~ Build from file ~~~~~~~~~~~~~~~~"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
ls -lhS $big_file
$time $exe build -k 20 -o graph $big_file
$time $exe stats graph --count-dummy

echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "~~~~~~~~~~~~ Build from file parallel ~~~~~~~~~~~"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
ls -lhS $big_file
$time $exe build -p $NUM_THREADS -k 20 -o graph $big_file
$time $exe stats graph -p $NUM_THREADS --count-dummy

echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "~~~~~~~~~~~~~~~~~ Build one big ~~~~~~~~~~~~~~~~~"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
$time $exe build -k 20 -o graph $DATA
$time $exe stats graph --count-dummy

echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "~~~~~~~~~~~~ Build one big parallel ~~~~~~~~~~~~~"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
ls -1S $DATA | $time $exe build -p $NUM_THREADS -k 20 -o graph
$time $exe stats graph -p $NUM_THREADS --count-dummy

echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "~~~~~~~~~~~~ Count k-mers parallel ~~~~~~~~~~~~~~"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
$time $(
for x in $DATA; do
    ./KMC/kmc -ci1 -t$NUM_THREADS -k20 -m5 -fa -b $x ${x}.kmc ./KMC
done
)

echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "~~~~~~~~~~~ Build from KMC parallel  ~~~~~~~~~~~~"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
ls -1S $DATA_PATH/*.kmc.kmc_pre | $time $exe build -p $NUM_THREADS -k 20 -o graph --kmc

echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "~~~~~~~~~~~ Build separate from KMC  ~~~~~~~~~~~~"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
$time $(
for x in $DATA_PATH/*.kmc.kmc_pre; do
    $exe build -k 20 -o ${x%.kmc.kmc_pre} --kmc $x
done
)

echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "~~~~~~~~~~~~~~~ Extract contigs  ~~~~~~~~~~~~~~~~"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
$time $(
for x in $DATA_PATH/*.dbg; do
    $exe assemble -o ${x%.dbg}.contigs $x
done
)

echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "~~~~~~~~~~~~~ Build from contigs  ~~~~~~~~~~~~~~~"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
ls -1S $DATA_PATH/*.contigs* | $time $exe build -k 20 -o graph

echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "~~~~~~~~~~~~~~~~ Merge graphs  ~~~~~~~~~~~~~~~~~~"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
ls -1S $DATA_PATH/*.dbg | $time $exe merge -o graph
