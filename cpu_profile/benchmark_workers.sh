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
WORKERS=(1 2 4 8 16 32 64 128 256)

# === Build ===
echo "[+] Building project..."
cd ..
mkdir -p "$SCRIPT_DIR/../$BUILD_DIR"
cd "$SCRIPT_DIR/../$BUILD_DIR"
cmake .. > /dev/null
make -j > /dev/null
cd cpp

# === CSV Header ===
echo "Workers,Runtime(s),IPC,CacheMisses,L1Misses,L2Misses,L3Misses,BackendStalls,MemoryStalls,BranchMisses,CyclesPerInstruction,MemoryBandwidth(MB/s)" > "$OUTPUT_CSV"

# === Benchmark Loop ===
for W in "${WORKERS[@]}"; do
    echo "[+] Running with CILK_NWORKERS=$W"

    N_OUT="$SCRIPT_DIR/nbody_out.txt"
    P_STAT="$SCRIPT_DIR/perf_stat.txt"

    # Remove temp files before run
    rm -f perf.data perf_script.data "$N_OUT" "$P_STAT"

    # === Run simulation silently ===
    CILK_NWORKERS=$W perf record -e cpu-cycles,instructions -c 3333333 -g -- ./$EXECUTABLE -n $BODY_COUNT > "$N_OUT" 2>/dev/null || continue

    # === Collect performance stats silently ===
    CILK_NWORKERS=$W perf stat -e \
        cpu-cycles,instructions,\
        cache-misses,cache-references,\
        L1-dcache-load-misses,L1-dcache-loads,\
        L2-dcache-load-misses,L2-dcache-loads,\
        LLC-load-misses,LLC-loads,\
        cpu/event=0xa3,umask=0x01,name=RESOURCE_STALLS.ANY/,\
        cpu/event=0x0d,umask=0x01,name=MEMORY_STALLS.ANY/,\
        branch-misses,branch-instructions,\
        mem-loads,mem-stores \
        -r 1 ./$EXECUTABLE -n $BODY_COUNT 1>/dev/null 2> "$P_STAT"

    # === Extract values ===
    RUNTIME=$(grep "Total time" "$N_OUT" | awk '{print $4}')
    INSTR=$(grep -i instructions "$P_STAT" | awk '{print $1}' | sed 's/,//g')
    CYCLES=$(grep -i cpu-cycles "$P_STAT" | awk '{print $1}' | sed 's/,//g')
    L1_MISS=$(grep -i L1-dcache-load-misses "$P_STAT" | awk '{print $1}' | sed 's/,//g')
    L2_MISS=$(grep -i L2-dcache-load-misses "$P_STAT" | awk '{print $1}' | sed 's/,//g')
    L3_MISS=$(grep -i LLC-load-misses "$P_STAT" | awk '{print $1}' | sed 's/,//g')
    STALLS=$(grep -i RESOURCE_STALLS.ANY "$P_STAT" | awk '{print $1}' | sed 's/,//g')
    MEM_STALLS=$(grep -i MEMORY_STALLS.ANY "$P_STAT" | awk '{print $1}' | sed 's/,//g')
    BRANCH_MISS=$(grep -i branch-misses "$P_STAT" | awk '{print $1}' | sed 's/,//g')
    MEM_LOADS=$(grep -i mem-loads "$P_STAT" | awk '{print $1}' | sed 's/,//g')
    MEM_STORES=$(grep -i mem-stores "$P_STAT" | awk '{print $1}' | sed 's/,//g')

    if [[ -n "$CYCLES" && "$CYCLES" != "0" ]]; then
        IPC=$(awk "BEGIN {printf \"%.3f\", $INSTR / $CYCLES}")
        CPI=$(awk "BEGIN {printf \"%.3f\", $CYCLES / $INSTR}")
    else
        IPC=""
        CPI=""
    fi

    # Calculate memory bandwidth (MB/s)
    if [[ -n "$RUNTIME" && -n "$MEM_LOADS" && -n "$MEM_STORES" ]]; then
        # Assuming 64 bytes per memory operation
        MEM_BW=$(awk "BEGIN {printf \"%.2f\", ($MEM_LOADS + $MEM_STORES) * 64 / $RUNTIME / 1000000}")
    else
        MEM_BW=""
    fi

    if [[ -n "$RUNTIME" && -n "$IPC" && -n "$L1_MISS" && -n "$L2_MISS" && -n "$L3_MISS" ]]; then
        echo "$W,$RUNTIME,$IPC,$L3_MISS,$L1_MISS,$L2_MISS,$L3_MISS,$STALLS,$MEM_STALLS,$BRANCH_MISS,$CPI,$MEM_BW" >> "$OUTPUT_CSV"
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