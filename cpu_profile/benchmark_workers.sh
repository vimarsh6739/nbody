#!/bin/bash

set -e

# === Config ===
SCRIPT_DIR=$(dirname "$(realpath "$0")")
BUILD_DIR=build
EXECUTABLE=cpu_nbody
FLAMEGRAPH_DIR=FlameGraph
OUTPUT_CSV="$SCRIPT_DIR/benchmark_results.csv"
PLOT_SCRIPT="$SCRIPT_DIR/plot_benchmarks.py"
BODY_COUNT=10000
WORKERS=(1 2 4 8 16 32 64 128)

# === Build ===
echo "[+] Building project..."
cd ..
mkdir -p "$SCRIPT_DIR/../$BUILD_DIR"
cd "$SCRIPT_DIR/../$BUILD_DIR"
cmake .. > /dev/null
make -j > /dev/null
cd cpp

# === CSV Header ===
echo "Workers,Runtime(s),IPC,CacheMisses,BackendStalls" > "$OUTPUT_CSV"

# === Benchmark Loop ===
for W in "${WORKERS[@]}"; do
    echo "[+] Running with CILK_NWORKERS=$W"

    N_OUT="$SCRIPT_DIR/nbody_out.txt"
    P_STAT="$SCRIPT_DIR/perf_stat.txt"

    # Remove temp files before run
    rm -f perf.data perf_script.data "$N_OUT" "$P_STAT"

    # === Run the simulation and collect perf record ===
    CILK_NWORKERS=$W perf record -e cpu-cycles,instructions -c 3333333 -g -- ./$EXECUTABLE -n $BODY_COUNT > "$N_OUT" 2>&1 || continue

    # === Collect performance stats ===
    CILK_NWORKERS=$W perf stat -e cpu-cycles,instructions,cache-misses,cpu/event=0xa3,umask=0x01,name=RESOURCE_STALLS.ANY/ -r 1 \
        ./$EXECUTABLE -n $BODY_COUNT 2> "$P_STAT"

    # === Extract values ===
    RUNTIME=$(grep "Total time" "$N_OUT" | awk '{print $4}')
    INSTR=$(grep -i instructions "$P_STAT" | awk '{print $1}' | sed 's/,//g')
    CYCLES=$(grep -i cpu-cycles "$P_STAT" | awk '{print $1}' | sed 's/,//g')
    MISS=$(grep -i cache-misses "$P_STAT" | awk '{print $1}' | sed 's/,//g')
    STALLS=$(grep -i RESOURCE_STALLS.ANY "$P_STAT" | awk '{print $1}' | sed 's/,//g')

    if [[ -n "$CYCLES" && "$CYCLES" != "0" ]]; then
        IPC=$(awk "BEGIN {printf \"%.3f\", $INSTR / $CYCLES}")
    else
        IPC=""
    fi

    if [[ -n "$RUNTIME" && -n "$IPC" && -n "$MISS" && -n "$STALLS" ]]; then
        echo "$W,$RUNTIME,$IPC,$MISS,$STALLS" >> "$OUTPUT_CSV"
    else
        echo "[!] Skipped logging for workers=$W due to missing data"
    fi

    # Clean up temp files after each run
    rm -f perf.data perf_script.data "$N_OUT" "$P_STAT"

done

# === Run Python Plot ===
echo "[+] Generating plots from $OUTPUT_CSV ..."
python3 "$PLOT_SCRIPT" "$OUTPUT_CSV"

echo "[✔] Benchmarking complete. Results saved to:"
echo " - CSV: $OUTPUT_CSV"
echo " - Plots: $SCRIPT_DIR/*.png"