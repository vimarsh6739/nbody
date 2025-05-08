import pandas as pd
import matplotlib.pyplot as plt
import sys
import os

# === Load CSV ===
if len(sys.argv) > 1:
    csv_path = sys.argv[1]
else:
    csv_path = os.path.join(os.path.dirname(__file__), 'benchmark_results.csv')

df = pd.read_csv(csv_path)

# === Ensure Workers column is sorted and treated as categorical
df = df.sort_values(by='Workers')
df['Workers'] = df['Workers'].astype(int)

# === Plot Config ===
metrics = ['Runtime(s)', 'IPC', 'CacheMisses', 'BackendStalls']
ylabels = ['Runtime (s)', 'IPC', 'Cache Misses', 'Backend Stalls']
xticks = [1, 2, 4, 8, 16, 32, 64, 128]

# === Plot each metric ===
for metric, ylabel in zip(metrics, ylabels):
    plt.figure()
    plt.plot(df['Workers'], df[metric], marker='o')
    plt.xlabel('CILK_NWORKERS')
    plt.ylabel(ylabel)
    plt.title(f'{ylabel} vs Workers')
    plt.xticks(xticks)  # Set specific x-axis ticks
    plt.grid(True)
    output_path = os.path.join(os.path.dirname(csv_path), f"{metric}.png")
    plt.savefig(output_path)
    print(f"[✔] Saved plot: {output_path}")