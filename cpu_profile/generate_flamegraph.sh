#!/bin/bash

set -e

# === Config ===
SCRIPT_DIR=$(dirname "$(realpath "$0")")
BUILD_DIR=build
EXECUTABLE=cpu_nbody
FLAMEGRAPH_DIR=FlameGraph
CILK_WORKERS=${1:-4}
NUM_BODIES=2000

# === Parse arguments ===
while [[ $# -gt 0 ]]; do
    case "$1" in
        -n|--nbody)
            NUM_BODIES="$2"
            shift 2
            ;;
        *)
            CILK_WORKERS="$1"
            shift
            ;;
    esac
done

# === Build ===
cd "$SCRIPT_DIR/.."
echo "[+] Building project..."
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"
cmake .. > /dev/null
make -j > /dev/null

# === Run perf ===
cd cpp
echo "[+] Running $EXECUTABLE with CILK_NWORKERS=$CILK_WORKERS and $NUM_BODIES bodies..."
rm -f perf.data perf_script.data folded_cycles folded_insts out.perf out.folded

# === Standard FlameGraph ===
CILK_NWORKERS=$CILK_WORKERS perf record -e cpu-cycles,instructions -c 3333333 -g ./$EXECUTABLE -n $NUM_BODIES
perf script > perf_script.data

echo "[+] Generating standard flamegraph..."
cp perf_script.data out.perf
"$SCRIPT_DIR/../$FLAMEGRAPH_DIR/stackcollapse-perf.pl" out.perf > out.folded
"$SCRIPT_DIR/../$FLAMEGRAPH_DIR/flamegraph.pl" out.folded > "$SCRIPT_DIR/flamegraph.svg"

# === IPC FlameGraph ===
echo "[+] Generating IPC flamegraph..."
"$SCRIPT_DIR/../$FLAMEGRAPH_DIR/stackcollapse-perf.pl" --event-filter=cpu-cycles perf_script.data > folded_cycles
"$SCRIPT_DIR/../$FLAMEGRAPH_DIR/stackcollapse-perf.pl" --event-filter=instructions perf_script.data > folded_insts

"$SCRIPT_DIR/../$FLAMEGRAPH_DIR/difffolded.pl" folded_insts folded_cycles \
| "$SCRIPT_DIR/../$FLAMEGRAPH_DIR/flamegraph.pl" \
    --countname=cycles \
    --title="IPC Flamegraph (red = IPC &lt; 1)" \
    --subtitle="Differential: instructions vs cycles" > "$SCRIPT_DIR/ipc_flamegraph.svg"

# === Backend Stall FlameGraph (Intel) ===
echo "[+] Generating backend-stall flamegraph..."
CILK_NWORKERS=$CILK_WORKERS perf record -F 999 -g \
    -e cycles:pp -e cpu/event=0xa3,umask=0x01,name=RESOURCE_STALLS.ANY/ \
    -- ./$EXECUTABLE -n $NUM_BODIES

perf script | "$SCRIPT_DIR/../$FLAMEGRAPH_DIR/stackcollapse-perf.pl" --all | \
    "$SCRIPT_DIR/../$FLAMEGRAPH_DIR/flamegraph.pl" --title "Backend stall cycles" \
                                        --color=hot \
                                        --countname="cycles" \
                                        > "$SCRIPT_DIR/backend-stall.svg"

# === Cache-miss FlameGraph ===
echo "[+] Generating cache-miss flamegraph..."
perf record -F 999 -g \
    -e cache-misses,cache-references \
    -- ./$EXECUTABLE -n $NUM_BODIES

perf script | "$SCRIPT_DIR/../$FLAMEGRAPH_DIR/stackcollapse-perf.pl" --all | \
    "$SCRIPT_DIR/../$FLAMEGRAPH_DIR/flamegraph.pl" \
        --title="LLC&amp;DRAM load misses" \
        --color=mem --countname="misses" \
        > "$SCRIPT_DIR/cache-miss.svg"

# === Cleanup temporary files ===
rm -f perf.data perf.data.old perf_script.data out.perf out.folded folded_cycles folded_insts

echo
echo "[✔] Flamegraphs generated!"
echo " - Standard:      $SCRIPT_DIR/flamegraph.svg"
echo " - IPC:           $SCRIPT_DIR/ipc_flamegraph.svg"
echo " - Backend Stall: $SCRIPT_DIR/backend-stall.svg"
echo " - Cache Miss:    $SCRIPT_DIR/cache-miss.svg"