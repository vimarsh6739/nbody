#!/bin/bash

set -e

# === Config ===
BUILD_DIR=build
EXECUTABLE=cpu_nbody
CILK_WORKERS=${1:-4}
FLAMEGRAPH_DIR=FlameGraph

# === Check FlameGraph Tools ===
if [ ! -d "$FLAMEGRAPH_DIR" ]; then
    echo "[+] Cloning FlameGraph..."
    git clone https://github.com/brendangregg/FlameGraph.git
fi

# === Build ===
echo "[+] Building project..."
mkdir -p $BUILD_DIR
cd $BUILD_DIR
cmake ..
make -j

# === Run perf ===
cd src
echo "[+] Running $EXECUTABLE with CILK_NWORKERS=$CILK_WORKERS..."
rm -f perf.data perf_script.data flamegraph.svg ipc_flamegraph.svg folded_cycles folded_insts out.perf out.folded

CILK_NWORKERS=$CILK_WORKERS \
perf record -e cpu-cycles,instructions -c 3333333 -g ./$EXECUTABLE

perf script > perf_script.data

# === Generate standard FlameGraph ===
echo "[+] Generating standard flamegraph..."
cp perf_script.data out.perf
../../$FLAMEGRAPH_DIR/stackcollapse-perf.pl out.perf > out.folded
../../$FLAMEGRAPH_DIR/flamegraph.pl out.folded > flamegraph.svg

# === Generate IPC FlameGraph ===
echo "[+] Generating IPC flamegraph..."
../../$FLAMEGRAPH_DIR/stackcollapse-perf.pl --event-filter=cpu-cycles perf_script.data > folded_cycles
../../$FLAMEGRAPH_DIR/stackcollapse-perf.pl --event-filter=instructions perf_script.data > folded_insts

../../$FLAMEGRAPH_DIR/difffolded.pl folded_insts folded_cycles \
| ../../$FLAMEGRAPH_DIR/flamegraph.pl \
    --countname=cycles \
    --title="IPC Flamegraph (red = IPC &lt; 1)" \
    --subtitle="Differential: instructions vs cycles" > ipc_flamegraph.svg

echo
echo "[✔] Flamegraphs generated!"
echo " - Standard: $(pwd)/flamegraph.svg"
echo " - IPC:      $(pwd)/ipc_flamegraph.svg"