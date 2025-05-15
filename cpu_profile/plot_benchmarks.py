#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os

def plot_metrics(df, output_dir):
    # Set style
    sns.set_style("whitegrid")
    plt.rcParams['figure.figsize'] = [12, 8]
    
    # 1. Strong Scaling Plot
    plt.figure()
    plt.plot(df['Workers'], df['Runtime(s)'], 'b-o', label='Actual Runtime')
    plt.plot(df['Workers'], df['Runtime(s)'][0] / df['Workers'], 'r--', label='Ideal Scaling')
    plt.xlabel('Number of Workers')
    plt.ylabel('Runtime (seconds)')
    plt.title('Strong Scaling Performance')
    plt.legend()
    plt.savefig(os.path.join(output_dir, 'strong_scaling.png'))
    plt.close()

    # 2. IPC and CPI
    plt.figure()
    ax1 = plt.gca()
    ax2 = ax1.twinx()
    ax1.plot(df['Workers'], df['IPC'], 'b-o', label='IPC')
    ax2.plot(df['Workers'], df['CyclesPerInstruction'], 'r-o', label='CPI')
    ax1.set_xlabel('Number of Workers')
    ax1.set_ylabel('Instructions Per Cycle (IPC)', color='b')
    ax2.set_ylabel('Cycles Per Instruction (CPI)', color='r')
    plt.title('Instruction Pipeline Efficiency')
    plt.savefig(os.path.join(output_dir, 'pipeline_efficiency.png'))
    plt.close()

    # 3. Cache Hierarchy Performance
    plt.figure()
    plt.plot(df['Workers'], df['L1Misses'], 'b-o', label='L1 Cache')
    plt.plot(df['Workers'], df['L2Misses'], 'r-o', label='L2 Cache')
    plt.plot(df['Workers'], df['L3Misses'], 'g-o', label='L3 Cache')
    plt.xlabel('Number of Workers')
    plt.ylabel('Cache Misses')
    plt.title('Cache Miss Rates Across Hierarchy')
    plt.legend()
    plt.savefig(os.path.join(output_dir, 'cache_performance.png'))
    plt.close()

    # 4. Memory and Branch Performance
    plt.figure()
    ax1 = plt.gca()
    ax2 = ax1.twinx()
    ax1.plot(df['Workers'], df['MemoryBandwidth(MB/s)'], 'b-o', label='Memory Bandwidth')
    ax2.plot(df['Workers'], df['BranchMisses'], 'r-o', label='Branch Misses')
    ax1.set_xlabel('Number of Workers')
    ax1.set_ylabel('Memory Bandwidth (MB/s)', color='b')
    ax2.set_ylabel('Branch Misses', color='r')
    plt.title('Memory and Branch Prediction Performance')
    plt.savefig(os.path.join(output_dir, 'memory_branch_performance.png'))
    plt.close()

    # 5. Stalls Analysis
    plt.figure()
    plt.plot(df['Workers'], df['BackendStalls'], 'b-o', label='Backend Stalls')
    plt.plot(df['Workers'], df['MemoryStalls'], 'r-o', label='Memory Stalls')
    plt.xlabel('Number of Workers')
    plt.ylabel('Number of Stalls')
    plt.title('Pipeline Stalls Analysis')
    plt.legend()
    plt.savefig(os.path.join(output_dir, 'stalls_analysis.png'))
    plt.close()

def main():
    if len(sys.argv) != 2:
        print("Usage: python3 plot_benchmarks.py <csv_file>")
        sys.exit(1)

    csv_file = sys.argv[1]
    output_dir = os.path.dirname(csv_file)
    
    # Read data
    df = pd.read_csv(csv_file)
    
    # Generate plots
    plot_metrics(df, output_dir)
    print(f"Plots generated in {output_dir}")

if __name__ == "__main__":
    main()